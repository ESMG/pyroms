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
from pyroms import _remapping

import matplotlib.pyplot as plt

def remapping(varname, srcfile, wts_files, srcgrd, dstgrd, \
              rotate_uv=False, trange=None, irange=None, jrange=None, \
              dstdir='./' ,zlevel=None, dmax=0, cdepth=0, kk=0, \
              uvar='u', vvar='v', rotate_part=False):
    '''
    A remapping function to go from a ROMS grid to another ROMS grid.
    If the u/v variables need to be rotated, it must be called for each
    u/v pair (such as u/v, uice/vice).
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

    # build intermediate zgrid
    if zlevel is None:
        zlevel = np.array([-7500.,-7000.,-6500.,-6000.,-5500.,-5000.,\
                   -4500.,-4000.,-3500.,-3000.,-2500.,-2000.,-1750.,\
                   -1500.,-1250.,-1000.,-900.,-800.,-700.,-600.,-500.,\
                   -400.,-300.,-250.,-200.,-175.,-150.,-125.,-100.,-90.,\
                   -80.,-70.,-60.,-50.,-45.,-40.,-35.,-30.,-25.,-20.,-17.5,\
                   -15.,-12.5,-10.,-7.5,-5.,-2.5,0.])
    else:
        zlevel = np.sort(-abs(zlevel))
    nzlevel = len(zlevel)
    src_zcoord = pyroms.vgrid.z_coordinate(srcgrd.vgrid.h, zlevel, nzlevel)
    dst_zcoord = pyroms.vgrid.z_coordinate(dstgrd.vgrid.h, zlevel, nzlevel)
    srcgrdz = pyroms.grid.ROMS_Grid(srcgrd.name+'_Z', srcgrd.hgrid, src_zcoord)
    dstgrdz = pyroms.grid.ROMS_Grid(dstgrd.name+'_Z', dstgrd.hgrid, dst_zcoord)

    # varname argument
    if type(varname).__name__ == 'list':
        nvar = len(varname)
    elif type(varname).__name__ == 'str':
        varname = [varname]
        nvar = len(varname)
    else:
        raise ValueError('varname must be a str or a list of str')

    # if we're working on u and v, we'll compute ubar,vbar afterwards
    compute_ubar = False
    if (varname.__contains__('u') == 1 and varname.__contains__('v') == 1) or \
       (varname.__contains__('u_eastward') == 1 and varname.__contains__('v_northward') == 1):
        compute_ubar = True
        print('ubar/vbar to be computed from u/v')
        if varname.__contains__('ubar'):
            varname.remove('ubar')
            nvar = nvar-1
        if varname.__contains__('vbar'):
            varname.remove('vbar')
            nvar = nvar-1

    # if rotate_uv=True, check that u and v are in varname
    if rotate_uv is True:
        if varname.__contains__(uvar) == 0 or varname.__contains__(vvar) == 0:
            raise Warning('varname must include uvar and vvar in order to' \
                   + ' rotate the velocity field')
        else:
            varname.remove(uvar)
            varname.remove(vvar)
            nvar = nvar-2

    # srcfile argument
    print('files', srcfile)
    if type(srcfile).__name__ == 'list':
        nfile = len(srcfile)
    elif type(srcfile).__name__ == 'str':
        srcfile = sorted(glob.glob(srcfile))
        nfile = len(srcfile)
    else:
        raise ValueError('src_srcfile must be a str or a list of str')
    print('number of files', nfile, srcfile)

    # get wts_file
    if type(wts_files).__name__ == 'str':
        wts_files = sorted(glob.glob(wts_files))

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
                     + dstgrd.name + '.nc'
            if os.path.exists(dstfile) is False:
                print('Creating destination file', dstfile)
                pyroms_toolbox.nc_create_roms_file(dstfile, dstgrd, ocean_time)

            # open destination file
            nc = netCDF.Dataset(dstfile, 'a', format='NETCDF3_64BIT')

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

                # determine variable dimension
                ndim = len(src_var.dimensions)-1

                # get spval
                try:
                    spval = src_var._FillValue
                except:
#                    raise Warning, 'Did not find a _FillValue attribute.'
                    print('Warning, Did not find a _FillValue attribute.')
                    spval = 1.e37

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
                if src_var.dimensions[2].find('_u') != -1:
                    Cpos='u'
                if src_var.dimensions[2].find('_v') != -1:
                    Cpos='v'
                if src_var.dimensions[1].find('_w') != -1:
                    Cpos='w'

                print('Arakawa C-grid position is', Cpos)

                # create variable in _destination file
                if nctidx == 0:
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

                # get the right remap weights file
                for s in range(len(wts_files)):
                    if wts_files[s].__contains__(Cpos+'_to_'+Cpos+'.nc'):
                        wts_file = wts_files[s]
                        break
                    else:
                        if s == len(wts_files) - 1:
                            raise ValueError('Did not find the appropriate remap weights file')

                if ndim == 3:
                    # vertical interpolation from sigma to standard z level
                    print('vertical interpolation from sigma to standard z level')
                    src_varz = pyroms.remapping.roms2z( \
                                 src_var[nt,:,jjrange[0]:jjrange[1],iirange[0]:iirange[1]], \
                                 srcgrd, srcgrdz, Cpos=Cpos, spval=spval, \
                                 irange=iirange, jrange=jjrange)

                    # flood the grid
                    print('flood the grid')
                    src_varz = pyroms.remapping.flood(src_varz, srcgrdz, Cpos=Cpos, \
                                      irange=iirange, jrange=jjrange, spval=spval, \
                                      dmax=dmax, cdepth=cdepth, kk=kk)

                else:
                    src_varz = src_var[nt,jjrange[0]:jjrange[1],iirange[0]:iirange[1]]

                # horizontal interpolation using scrip weights
                print('horizontal interpolation using scrip weights')
                dst_varz = pyroms.remapping.remap(src_varz, wts_file, \
                                                  spval=spval)


                if ndim == 3:
                    # vertical interpolation from standard z level to sigma
                    print('vertical interpolation from standard z level to sigma')
                    dst_var = pyroms.remapping.z2roms(dst_varz, dstgrdz, dstgrd, \
                                     Cpos=Cpos, spval=spval, flood=False)
                else:
                    dst_var = dst_varz

                if varname[nv] == 'u':
                    dst_u = dst_var
                if varname[nv] == 'v':
                    dst_v = dst_var

                # write data in destination file
                print('write data in destination file')
                nc.variables[varname[nv]][nctidx] = dst_var


            # rotate the velocity field if requested
            if rotate_uv is True:
                print(' ')
                print('remapping and rotating', uvar, 'and', vvar, 'from', \
                      srcgrd.name, 'to', dstgrd.name)

                # get source data
                src_u = pyroms.utility.get_nc_var(uvar, srcfile[nf])
                src_v = pyroms.utility.get_nc_var(vvar, srcfile[nf])

                # get spval
                try:
                    spval = src_v._FillValue
                except:
                    raise Warning('Did not find a _FillValue attribute.')

                if rotate_part:
                    ndim = len(src_u.dimensions)-1
                    ind = uvar.find('_eastward')
                    uvar_out = uvar[0:ind]
                    print("Warning: renaming uvar to", uvar_out)
#                   print("uvar dims:", src_u.dimensions)
                    ind = vvar.find('_northward')
                    vvar_out = vvar[0:ind]
                    print("Warning: renaming vvar to", vvar_out)
#                   print("vvar dims:", src_v.dimensions)
                    if ndim == 3:
                        dimens_u = ['ocean_time', 's_rho', 'eta_u', 'xi_u']
                        dimens_v = ['ocean_time', 's_rho', 'eta_v', 'xi_v']
                    else:
                        dimens_u = ['ocean_time', 'eta_u', 'xi_u']
                        dimens_v = ['ocean_time', 'eta_v', 'xi_v']

                else:
                    dimens_u = [i for i in src_u.dimensions]
                    dimens_v = [i for i in src_v.dimensions]
                    uvar_out = uvar
                    vvar_out = vvar

                # create variable in destination file
                if nctidx == 0:
                    print('Creating variable '+uvar_out)
                    nc.createVariable(uvar_out, 'f8', dimens_u, fill_value=spval)
                    nc.variables[uvar_out].long_name = src_u.long_name
                    nc.variables[uvar_out].units = src_u.units
                    nc.variables[uvar_out].time = src_u.time
                    nc.variables[uvar_out].coordinates = \
                           str(dimens_u.reverse())
                    nc.variables[uvar_out].field = src_u.field
                    print('Creating variable '+vvar_out)
                    nc.createVariable(vvar_out, 'f8', dimens_v, fill_value=spval)
                    nc.variables[vvar_out].long_name = src_v.long_name
                    nc.variables[vvar_out].units = src_v.units
                    nc.variables[vvar_out].time = src_v.time
                    nc.variables[vvar_out].coordinates = \
                           str(dimens_v.reverse())
                    nc.variables[vvar_out].field = src_v.field

                # get the right remap weights file
                if rotate_part:
                    for s in range(len(wts_files)):
                        if wts_files[s].__contains__('rho_to_rho.nc'):
                            wts_file_u = wts_files[s]
                            wts_file_v = wts_files[s]
                    Cpos_u = 'rho'
                    Cpos_v = 'rho'
                else:
                    for s in range(len(wts_files)):
                        if wts_files[s].__contains__('u_to_rho.nc'):
                            wts_file_u = wts_files[s]
                        if wts_files[s].__contains__('v_to_rho.nc'):
                            wts_file_v = wts_files[s]
                    Cpos_u = 'u'
                    Cpos_v = 'v'

                # get the right ranges
                if rotate_part:
                    # irange
                    if irange is None:
                        iirange = (0,src_u.shape[-1])
                    else:
                        iirange = irange
                    # jrange
                    if jrange is None:
                        jjrange = (0,src_u.shape[-2])
                    else:
                        jjrange = jrange
                else:
                    # irange
                    if irange is None:
                        iirange = (0,src_u.shape[-1])
                    else:
                        iirange = (irange[0], irange[1]-1)
                    # jrange
                    if jrange is None:
                        jjrange = (0,src_u.shape[-2])
                    else:
                        jjrange = jrange

                # vertical interpolation from sigma to standard z level

                ndim = len(src_v.dimensions)-1
                if ndim == 3:
                    print('vertical interpolation from sigma to standard z level')
                    src_uz = pyroms.remapping.roms2z( \
                            src_u[nt,:,jjrange[0]:jjrange[1],iirange[0]:iirange[1]], \
                            srcgrd, srcgrdz, Cpos=Cpos_u, spval=spval, \
                            irange=iirange, jrange=jjrange)
                    # flood the grid
                    print('flood the u grid')
                    src_uz = pyroms.remapping.flood(src_uz, srcgrdz, Cpos=Cpos_u, \
                                  irange=iirange, jrange=jjrange, \
                                  spval=spval, dmax=dmax, cdepth=cdepth, kk=kk)
                else:
                    src_uz = src_u[nt,jjrange[0]:jjrange[1],iirange[0]:iirange[1]]
                    src_uz = pyroms.remapping.flood2d(src_uz, srcgrdz, Cpos=Cpos_u, \
                                      irange=iirange, jrange=jjrange, spval=spval, \
                                      dmax=dmax)

                # get the right ranges
                if rotate_part:
                    # irange
                    if irange is None:
                        iirange = (0,src_v.shape[-1])
                    else:
                        iirange = irange
                    # jrange
                    if jrange is None:
                        jjrange = (0,src_v.shape[-2])
                    else:
                        jjrange = jrange
                else:
                    # irange
                    if irange is None:
                        iirange = (0,src_v.shape[-1])
                    else:
                        iirange = irange
                    # jrange
                    if jrange is None:
                        jjrange = (0,src_v.shape[-2])
                    else:
                        jjrange = (jrange[0], jrange[1]-1)

                if ndim == 3:
                    src_vz = pyroms.remapping.roms2z( \
                            src_v[nt,:,jjrange[0]:jjrange[1],iirange[0]:iirange[1]], \
                            srcgrd, srcgrdz, Cpos=Cpos_v, spval=spval, \
                            irange=iirange, jrange=jjrange)

                    # flood the grid
                    print('flood the v grid')
                    src_vz = pyroms.remapping.flood(src_vz, srcgrdz, Cpos=Cpos_v, \
                                  irange=iirange, jrange=jjrange, \
                                  spval=spval, dmax=dmax, cdepth=cdepth, kk=kk)
                else:
                    src_vz = src_v[nt,jjrange[0]:jjrange[1],iirange[0]:iirange[1]]
                    src_vz = pyroms.remapping.flood2d(src_vz, srcgrdz, Cpos=Cpos_v, \
                                      irange=iirange, jrange=jjrange, spval=spval, \
                                      dmax=dmax)

                # horizontal interpolation using scrip weights
                print('horizontal interpolation using scrip weights')
                dst_uz = pyroms.remapping.remap(src_uz, wts_file_u, \
                                                  spval=spval)
                dst_vz = pyroms.remapping.remap(src_vz, wts_file_v, \
                                                  spval=spval)

                if ndim == 3:
                    # vertical interpolation from standard z level to sigma
                    print('vertical interpolation from standard z level to sigma')
                    dst_u = pyroms.remapping.z2roms(dst_uz, dstgrdz, dstgrd, \
                                 Cpos='rho', spval=spval, flood=False)
                    dst_v = pyroms.remapping.z2roms(dst_vz, dstgrdz, dstgrd, \
                                 Cpos='rho', spval=spval, flood=False)
                else:
                    dst_u = dst_uz
                    dst_v = dst_vz

                # rotate u,v fields
                if rotate_part:
                    src_angle = np.zeros(dstgrd.hgrid.angle_rho.shape)
                else:
                    for s in range(len(wts_files)):
                        if wts_files[s].__contains__('rho_to_rho.nc'):
                            wts_file = wts_files[s]
                    src_ang = srcgrd.hgrid.angle_rho[jjrange[0]:jjrange[1],iirange[0]:iirange[1]]
                    src_angle = pyroms.remapping.remap(src_ang, wts_file)

                dst_angle = dstgrd.hgrid.angle_rho
                angle = dst_angle - src_angle
                if ndim == 3:
                    angle = np.tile(angle, (dstgrd.vgrid.N, 1, 1))

                U = dst_u + dst_v*1j
                eitheta = np.exp(-1j*angle)
                U = U * eitheta

                dst_u = np.real(U)
                dst_v = np.imag(U)

                # spval
                idxu = np.where(dstgrd.hgrid.mask_u == 0)
                idxv = np.where(dstgrd.hgrid.mask_v == 0)

                # move back to u,v points
                if ndim == 3:
                    dst_u = 0.5 * (dst_u[:,:,:-1] + dst_u[:,:,1:])
                    dst_v = 0.5 * (dst_v[:,:-1,:] + dst_v[:,1:,:])
                    for n in range(dstgrd.vgrid.N):
                        dst_u[n,idxu[0], idxu[1]] = spval
                        dst_v[n,idxv[0], idxv[1]] = spval
                else:
                    dst_u = 0.5 * (dst_u[:,:-1] + dst_u[:,1:])
                    dst_v = 0.5 * (dst_v[:-1,:] + dst_v[1:,:])
                    dst_u[idxu[0], idxu[1]] = spval
                    dst_v[idxv[0], idxv[1]] = spval

                # write data in destination file
                print('write data in destination file')
                nc.variables[uvar_out][nctidx] = dst_u
                nc.variables[vvar_out][nctidx] = dst_v

            if compute_ubar:
                if nctidx == 0:
                    print('Creating variable ubar')
                    nc.createVariable('ubar', 'f8', \
                         ('ocean_time', 'eta_u', 'xi_u'), fill_value=spval)
                    nc.variables['ubar'].long_name = '2D u-momentum component'
                    nc.variables['ubar'].units = 'meter second-1'
                    nc.variables['ubar'].time = 'ocean_time'
                    nc.variables['ubar'].coordinates = 'xi_u eta_u ocean_time'
                    nc.variables['ubar'].field = 'ubar-velocity,, scalar, series'
                    print('Creating variable vbar')
                    nc.createVariable('vbar', 'f8', \
                         ('ocean_time', 'eta_v', 'xi_v'), fill_value=spval)
                    nc.variables['vbar'].long_name = '2D v-momentum component'
                    nc.variables['vbar'].units = 'meter second-1'
                    nc.variables['vbar'].time = 'ocean_time'
                    nc.variables['vbar'].coordinates = 'xi_v eta_v ocean_time'
                    nc.variables['vbar'].field = 'vbar-velocity,, scalar, series'

                # compute depth average velocity ubar and vbar
                # get z at the right position
                z_u = 0.5 * (dstgrd.vgrid.z_w[0,:,:,:-1] + \
                        dstgrd.vgrid.z_w[0,:,:,1:])
                z_v = 0.5 * (dstgrd.vgrid.z_w[0,:,:-1,:] + \
                        dstgrd.vgrid.z_w[0,:,1:,:])

                dst_ubar = np.zeros((dst_u.shape[1], dst_u.shape[2]))
                dst_vbar = np.zeros((dst_v.shape[1], dst_v.shape[2]))

                for i in range(dst_ubar.shape[1]):
                    for j in range(dst_ubar.shape[0]):
                        dst_ubar[j,i] = (dst_u[:,j,i] * \
                                np.diff(z_u[:,j,i])).sum() / -z_u[0,j,i]

                for i in range(dst_vbar.shape[1]):
                    for j in range(dst_vbar.shape[0]):
                        dst_vbar[j,i] = (dst_v[:,j,i] * \
                                np.diff(z_v[:,j,i])).sum() / -z_v[0,j,i]

                # spval
                idxu = np.where(dstgrd.hgrid.mask_u == 0)
                idxv = np.where(dstgrd.hgrid.mask_v == 0)
                dst_ubar[idxu[0], idxu[1]] = spval
                dst_vbar[idxv[0], idxv[1]] = spval

                nc.variables['ubar'][nctidx] = dst_ubar
                nc.variables['vbar'][nctidx] = dst_vbar

            nctidx = nctidx + 1
            print('ADDING to nctidx ', nctidx)
            nc.sync()

    # close destination file
    nc.close()

    return
