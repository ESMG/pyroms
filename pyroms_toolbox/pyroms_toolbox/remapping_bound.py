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

import datetime

def remapping_bound(varname, srcfile, wts_files, srcgrd, dst_grd, \
              rotate_uv=False, trange=None, irange=None, jrange=None, \
              dstdir='./' ,zlevel=None, dmax=0, cdepth=0, kk=0, \
              uvar='u', vvar='v', rotate_part=False):
    '''
    A remapping function to extract boundary conditions from one ROMS grid
    to another. It will optionally rotating u and v variables, but needs
    to be called separately for each u/v pair (such as u/v, uice/vice).
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
    dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
    srcgrdz = pyroms.grid.ROMS_Grid(srcgrd.name+'_Z', srcgrd.hgrid, src_zcoord)
    dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid, dst_zcoord)

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

                # determine variable dimension
                ndim = len(src_var.dimensions)-1

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
                if src_var.dimensions[2].find('_u') != -1:
                    Cpos='u'
                    Lp = Lp-1
                    if irange is not None:
                        iirange = (irange[0], irange[1]-1)
                if src_var.dimensions[2].find('_v') != -1:
                    Cpos='v'
                    Mp = Mp-1
                    if jrange is not None:
                        jjrange = (jrange[0], jrange[1]-1)
                if src_var.dimensions[1].find('_w') != -1:
                    Cpos='w'

                print('Arakawa C-grid position is', Cpos)

                # create variable in _destination file
                if nctidx == 0:
                    for sid in sides:
                       varn = varname[nv]+str(sid)
                       dimens = [i for i in src_var.dimensions]
                       for dim in dimens:
                           if re.match(dimexcl[sid],dim):
                               dimens.remove(dim)
                       print('Creating variable', varn, dimens)
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

                print(datetime.datetime.now())
                # horizontal interpolation using scrip weights
                print('horizontal interpolation using scrip weights')
                dst_varz = pyroms.remapping.remap(src_varz, wts_file, \
                                                  spval=spval)

                if ndim == 3:
                    dst_var_north = pyroms.remapping.z2roms(dst_varz[:, \
                          Mp-1:Mp,0:Lp], dst_grdz, dst_grd, Cpos=Cpos, \
                          spval=spval, flood=False, irange=(0,Lp), \
                          jrange=(Mp-1,Mp))
                    dst_var_south = pyroms.remapping.z2roms(dst_varz[:, \
                          0:1, :], dst_grdz, dst_grd, Cpos=Cpos, \
                          spval=spval, flood=False, irange=(0,Lp), \
                          jrange=(0,1))
                    dst_var_east = pyroms.remapping.z2roms(dst_varz[:, \
                          :, Lp-1:Lp], dst_grdz, dst_grd, Cpos=Cpos, \
                          spval=spval, flood=False, irange=(Lp-1,Lp), \
                          jrange=(0,Mp))
                    dst_var_west = pyroms.remapping.z2roms(dst_varz[:, \
                          :, 0:1], dst_grdz, dst_grd, Cpos=Cpos, \
                          spval=spval, flood=False, irange=(0,1), \
                          jrange=(0,Mp))
                    if varname[nv] == 'u':
                        dst_u_west = dst_var_west
                        dst_u_east = dst_var_east
                        dst_u_north = dst_var_north
                        dst_u_south = dst_var_south
                    if varname[nv] == 'v':
                        dst_v_west = dst_var_west
                        dst_v_east = dst_var_east
                        dst_v_north = dst_var_north
                        dst_v_south = dst_var_south

                else:
                    dst_var_north = dst_varz[-1, :]
                    dst_var_south = dst_varz[0, :]
                    dst_var_east = dst_varz[:, -1]
                    dst_var_west = dst_varz[:, 0]

#                print datetime.datetime.now()

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
            if rotate_uv is True:
                print(' ')
                print('remapping and rotating u and v from', srcgrd.name, \
                      'to', dst_grd.name)

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
                    ind = vvar.find('_northward')
                    vvar_out = vvar[0:ind]
                    print("Warning: renaming vvar to", vvar_out)
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
                    print('Creating boundary variables for '+uvar)
                    for sid in sides:
                       varn = uvar_out+str(sid)
                       print('Creating variable', varn)
                       dimens = list(dimens_u)
                       for dim in dimens:
                           if re.match(dimexcl[sid],dim):
                               dimens.remove(dim)
                       nc.createVariable(varn, 'f8', dimens, \
                         fill_value=spval)
                       nc.variables[varn].long_name = uvar_out + \
                           ' ' + long[sid] + ' boundary condition'
                       try:
                           nc.variables[varn].units = src_u.units
                       except:
                           print(varn+' has no units')
                       nc.variables[varn].time = src_u.time
                       nc.variables[varn].coordinates = \
                           str(dimens.reverse())
                       nc.variables[varn].field = src_u.field
                    print('Creating boundary variables for '+vvar)
                    for sid in sides:
                       varn = vvar_out+str(sid)
                       print('Creating variable', varn)
                       dimens = list(dimens_v)
                       for dim in dimens:
                           if re.match(dimexcl[sid],dim):
                               dimens.remove(dim)
                       nc.createVariable(varn, 'f8', dimens, \
                         fill_value=spval)
                       nc.variables[varn].long_name = vvar_out + \
                                ' ' + long[sid] + ' boundary condition'
                       try:
                           nc.variables[varn].units = src_v.units
                       except:
                           print(varn+' has no units')
                       nc.variables[varn].time = src_v.time
                       nc.variables[varn].coordinates = \
                           str(dimens.reverse())
                       nc.variables[varn].field = src_v.field

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
                Mp, Lp = dst_grd.hgrid.mask_rho.shape

                if ndim == 3:
                    # vertical interpolation from standard z level to sigma
                    print('vertical interpolation from standard z level to sigma')
                    dst_u_north = pyroms.remapping.z2roms(dst_uz[:, Mp-2:Mp, 0:Lp], \
                         dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                         flood=False, irange=(0,Lp), jrange=(Mp-2,Mp))
                    dst_u_south = pyroms.remapping.z2roms(dst_uz[:, 0:2, 0:Lp], \
                         dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                         flood=False, irange=(0,Lp), jrange=(0,2))
                    dst_u_east = pyroms.remapping.z2roms(dst_uz[:, 0:Mp, Lp-2:Lp], \
                         dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                         flood=False, irange=(Lp-2,Lp), jrange=(0,Mp))
                    dst_u_west = pyroms.remapping.z2roms(dst_uz[:, 0:Mp, 0:2], \
                         dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                         flood=False, irange=(0,2), jrange=(0,Mp))

                    dst_v_north = pyroms.remapping.z2roms(dst_vz[:, Mp-2:Mp, 0:Lp], \
                         dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                         flood=False, irange=(0,Lp), jrange=(Mp-2,Mp))
                    dst_v_south = pyroms.remapping.z2roms(dst_vz[:, 0:2, 0:Lp], \
                         dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                         flood=False, irange=(0,Lp), jrange=(0,2))
                    dst_v_east = pyroms.remapping.z2roms(dst_vz[:, 0:Mp, Lp-2:Lp], \
                         dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                         flood=False, irange=(Lp-2,Lp), jrange=(0,Mp))
                    dst_v_west = pyroms.remapping.z2roms(dst_vz[:, 0:Mp, 0:2], \
                      dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                      flood=False, irange=(0,2), jrange=(0,Mp))
                else:
                    dst_u_north = dst_uz[Mp-2:Mp, 0:Lp]
                    dst_u_south = dst_uz[0:2, 0:Lp]
                    dst_u_east = dst_uz[0:Mp, Lp-2:Lp]
                    dst_u_west = dst_uz[0:Mp, 0:2]
                    dst_v_north = dst_vz[Mp-2:Mp, 0:Lp]
                    dst_v_south = dst_vz[0:2, 0:Lp]
                    dst_v_east = dst_vz[0:Mp, Lp-2:Lp]
                    dst_v_west = dst_vz[0:Mp, 0:2]

                # rotate u,v fields
                if rotate_part:
                    src_angle = np.zeros(dst_grd.hgrid.angle_rho.shape)
                else:
                    for s in range(len(wts_files)):
                        if wts_files[s].__contains__('rho_to_rho.nc'):
                            wts_file = wts_files[s]
                    src_ang = srcgrd.hgrid.angle_rho[jjrange[0]:jjrange[1],iirange[0]:iirange[1]]
                    src_angle = pyroms.remapping.remap(src_ang, wts_file)

                dst_angle = dst_grd.hgrid.angle_rho
                angle = dst_angle - src_angle
                if ndim == 3:
                    angle = np.tile(angle, (dst_grd.vgrid.N, 1, 1))
                    U_north = dst_u_north + dst_v_north*1j
                    eitheta_north = np.exp(-1j*angle[:,Mp-2:Mp, 0:Lp])
                    U_south = dst_u_south + dst_v_south*1j
                    eitheta_south = np.exp(-1j*angle[:,0:2, 0:Lp])
                    U_east = dst_u_east + dst_v_east*1j
                    eitheta_east = np.exp(-1j*angle[:,0:Mp, Lp-2:Lp])
                    U_west = dst_u_west + dst_v_west*1j
                    eitheta_west = np.exp(-1j*angle[:,0:Mp, 0:2])
                else:
                    U_north = dst_u_north + dst_v_north*1j
                    eitheta_north = np.exp(-1j*angle[Mp-2:Mp, 0:Lp])
                    U_south = dst_u_south + dst_v_south*1j
                    eitheta_south = np.exp(-1j*angle[0:2, 0:Lp])
                    U_east = dst_u_east + dst_v_east*1j
                    eitheta_east = np.exp(-1j*angle[0:Mp, Lp-2:Lp])
                    U_west = dst_u_west + dst_v_west*1j
                    eitheta_west = np.exp(-1j*angle[0:Mp, 0:2])

                U_north = U_north * eitheta_north
                dst_u_north = np.real(U_north)
                dst_v_north = np.imag(U_north)

                U_south = U_south * eitheta_south
                dst_u_south = np.real(U_south)
                dst_v_south = np.imag(U_south)

                U_east = U_east * eitheta_east
                dst_u_east = np.real(U_east)
                dst_v_east = np.imag(U_east)

                U_west = U_west * eitheta_west
                dst_u_west = np.real(U_west)
                dst_v_east = np.imag(U_east)

                # move back to u,v points
                if ndim == 3:
                    dst_u_north = 0.5 * np.squeeze(dst_u_north[:,-1,:-1] + \
                            dst_u_north[:,-1,1:])
                    dst_v_north = 0.5 * np.squeeze(dst_v_north[:,:-1,:] + \
                            dst_v_north[:,1:,:])
                    dst_u_south = 0.5 * np.squeeze(dst_u_south[:,0,:-1] + \
                            dst_u_south[:,0,1:])
                    dst_v_south = 0.5 * np.squeeze(dst_v_south[:,:-1,:] + \
                            dst_v_south[:,1:,:])
                    dst_u_east = 0.5 * np.squeeze(dst_u_east[:,:,:-1] + \
                            dst_u_east[:,:,1:])
                    dst_v_east = 0.5 * np.squeeze(dst_v_east[:,:-1,-1] + \
                            dst_v_east[:,1:,-1])
                    dst_u_west = 0.5 * np.squeeze(dst_u_west[:,:,:-1] + \
                            dst_u_west[:,:,1:])
                    dst_v_west = 0.5 * np.squeeze(dst_v_west[:,:-1,0] + \
                            dst_v_west[:,1:,0])
                else:
                    dst_u_north = 0.5 * np.squeeze(dst_u_north[-1,:-1] + \
                            dst_u_north[-1,1:])
                    dst_v_north = 0.5 * np.squeeze(dst_v_north[:-1,:] + \
                            dst_v_north[1:,:])
                    dst_u_south = 0.5 * np.squeeze(dst_u_south[0,:-1] + \
                            dst_u_south[0,1:])
                    dst_v_south = 0.5 * np.squeeze(dst_v_south[:-1,:] + \
                            dst_v_south[1:,:])
                    dst_u_east = 0.5 * np.squeeze(dst_u_east[:,:-1] + \
                            dst_u_east[:,1:])
                    dst_v_east = 0.5 * np.squeeze(dst_v_east[:-1,-1] + \
                            dst_v_east[1:,-1])
                    dst_u_west = 0.5 * np.squeeze(dst_u_west[:,:-1] + \
                            dst_u_west[:,1:])
                    dst_v_west = 0.5 * np.squeeze(dst_v_west[:-1,0] + \
                            dst_v_west[1:,0])

                # spval
                idxu_north = np.where(dst_grd.hgrid.mask_u[-1,:] == 0)
                idxv_north = np.where(dst_grd.hgrid.mask_v[-1,:] == 0)
                idxu_south = np.where(dst_grd.hgrid.mask_u[0,:] == 0)
                idxv_south = np.where(dst_grd.hgrid.mask_v[0,:] == 0)
                idxu_east = np.where(dst_grd.hgrid.mask_u[:,-1] == 0)
                idxv_east = np.where(dst_grd.hgrid.mask_v[:,-1] == 0)
                idxu_west = np.where(dst_grd.hgrid.mask_u[:,0] == 0)
                idxv_west = np.where(dst_grd.hgrid.mask_v[:,0] == 0)
                if ndim == 3:
                    for n in range(dst_grd.vgrid.N):
                        dst_u_north[n, idxu_north[0]] = spval
                        dst_v_north[n, idxv_north[0]] = spval
                        dst_u_south[n, idxu_south[0]] = spval
                        dst_v_south[n, idxv_south[0]] = spval
                        dst_u_east[n, idxu_east[0]] = spval
                        dst_v_east[n, idxv_east[0]] = spval
                        dst_u_west[n, idxu_west[0]] = spval
                        dst_v_west[n, idxv_west[0]] = spval
                else:
                    dst_u_north[idxu_north[0]] = spval
                    dst_v_north[idxv_north[0]] = spval
                    dst_u_south[idxu_south[0]] = spval
                    dst_v_south[idxv_south[0]] = spval
                    dst_u_east[idxu_east[0]] = spval
                    dst_v_east[idxv_east[0]] = spval
                    dst_u_west[idxu_west[0]] = spval
                    dst_v_west[idxv_west[0]] = spval

                # write data in destination file
                print('write data in destination file')
                sid = '_west'
                varn = uvar_out+str(sid)
                nc.variables[varn][nctidx] = dst_u_west
                varn = vvar_out+str(sid)
                nc.variables[varn][nctidx] = dst_v_west

                sid = '_north'
                varn = uvar_out+str(sid)
                nc.variables[varn][nctidx] = dst_u_north
                varn = vvar_out+str(sid)
                nc.variables[varn][nctidx] = dst_v_north

                sid = '_east'
                varn = uvar_out+str(sid)
                nc.variables[varn][nctidx] = dst_u_east
                varn = vvar_out+str(sid)
                nc.variables[varn][nctidx] = dst_v_east

                sid = '_south'
                varn = uvar_out+str(sid)
                nc.variables[varn][nctidx] = dst_u_south
                varn = vvar_out+str(sid)
                nc.variables[varn][nctidx] = dst_v_south

            if compute_ubar:
                if nctidx == 0:
                    print('Creating variable ubar_north')
                    nc.createVariable('ubar_north', 'f8', \
                         ('ocean_time', 'xi_u'), fill_value=spval)
                    nc.variables['ubar_north'].long_name = \
                          '2D u-momentum north boundary condition'
                    nc.variables['ubar_north'].units = 'meter second-1'
                    nc.variables['ubar_north'].time = 'ocean_time'
                    nc.variables['ubar_north'].coordinates = 'xi_u ocean_time'
                    nc.variables['ubar_north'].field = 'ubar_north, scalar, series'
                    print('Creating variable vbar_north')
                    nc.createVariable('vbar_north', 'f8', \
                         ('ocean_time', 'xi_v'), fill_value=spval)
                    nc.variables['vbar_north'].long_name = \
                          '2D v-momentum north boundary condition'
                    nc.variables['vbar_north'].units = 'meter second-1'
                    nc.variables['vbar_north'].time = 'ocean_time'
                    nc.variables['vbar_north'].coordinates = 'xi_v ocean_time'
                    nc.variables['vbar_north'].field = 'vbar_north,, scalar, series'

                    print('Creating variable ubar_south')
                    nc.createVariable('ubar_south', 'f8', \
                         ('ocean_time', 'xi_u'), fill_value=spval)
                    nc.variables['ubar_south'].long_name = \
                          '2D u-momentum south boundary condition'
                    nc.variables['ubar_south'].units = 'meter second-1'
                    nc.variables['ubar_south'].time = 'ocean_time'
                    nc.variables['ubar_south'].coordinates = 'xi_u ocean_time'
                    nc.variables['ubar_south'].field = 'ubar_south, scalar, series'
                    print('Creating variable vbar_south')
                    nc.createVariable('vbar_south', 'f8', \
                         ('ocean_time', 'xi_v'), fill_value=spval)
                    nc.variables['vbar_south'].long_name = \
                          '2D v-momentum south boundary condition'
                    nc.variables['vbar_south'].units = 'meter second-1'
                    nc.variables['vbar_south'].time = 'ocean_time'
                    nc.variables['vbar_south'].coordinates = 'xi_v ocean_time'

                    print('Creating variable ubar_west')
                    nc.createVariable('ubar_west', 'f8', \
                         ('ocean_time', 'eta_u'), fill_value=spval)
                    nc.variables['ubar_west'].long_name = \
                          '2D u-momentum west boundary condition'
                    nc.variables['ubar_west'].units = 'meter second-1'
                    nc.variables['ubar_west'].time = 'ocean_time'
                    nc.variables['ubar_west'].coordinates = 'eta_u ocean_time'
                    nc.variables['ubar_west'].field = 'ubar_west, scalar, series'
                    print('Creating variable vbar_west')
                    nc.createVariable('vbar_west', 'f8', \
                         ('ocean_time', 'eta_v'), fill_value=spval)
                    nc.variables['vbar_west'].long_name = \
                          '2D v-momentum west boundary condition'
                    nc.variables['vbar_west'].units = 'meter second-1'
                    nc.variables['vbar_west'].time = 'ocean_time'
                    nc.variables['vbar_west'].coordinates = 'eta_v ocean_time'

                    print('Creating variable ubar_east')
                    nc.createVariable('ubar_east', 'f8', \
                         ('ocean_time', 'eta_u'), fill_value=spval)
                    nc.variables['ubar_east'].long_name = \
                          '2D u-momentum east boundary condition'
                    nc.variables['ubar_east'].units = 'meter second-1'
                    nc.variables['ubar_east'].time = 'ocean_time'
                    nc.variables['ubar_east'].coordinates = 'eta_u ocean_time'
                    nc.variables['ubar_east'].field = 'ubar_east, scalar, series'
                    print('Creating variable vbar_east')
                    nc.createVariable('vbar_east', 'f8', \
                         ('ocean_time', 'eta_v'), fill_value=spval)
                    nc.variables['vbar_east'].long_name = \
                          '2D v-momentum east boundary condition'
                    nc.variables['vbar_east'].units = 'meter second-1'
                    nc.variables['vbar_east'].time = 'ocean_time'
                    nc.variables['vbar_east'].coordinates = 'eta_v ocean_time'

                # compute depth average velocity ubar and vbar
                # get z at the right position
                print('Computing ubar/vbar from u/v')
                z_u_north = 0.5 * (dst_grd.vgrid.z_w[0,:,-1,:-1] +
                        dst_grd.vgrid.z_w[0,:,-1, 1:])
                z_v_north = 0.5 * (dst_grd.vgrid.z_w[0,:,-1,:] +
                        dst_grd.vgrid.z_w[0,:,-2,:])
                z_u_south = 0.5 * (dst_grd.vgrid.z_w[0,:,0,:-1] +
                        dst_grd.vgrid.z_w[0,:,0,1:])
                z_v_south = 0.5 * (dst_grd.vgrid.z_w[0,:,0,:] +
                        dst_grd.vgrid.z_w[0,:,1,:])
                z_u_east = 0.5 * (dst_grd.vgrid.z_w[0,:,:,-1] +
                        dst_grd.vgrid.z_w[0,:,:,-2])
                z_v_east = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,-1] +
                        dst_grd.vgrid.z_w[0,:,1:,-1])
                z_u_west = 0.5 * (dst_grd.vgrid.z_w[0,:,:,0] +
                        dst_grd.vgrid.z_w[0,:,:,1])
                z_v_west = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,0] +
                        dst_grd.vgrid.z_w[0,:,1:,0])
                if not rotate_uv:
                    dst_u_north = np.squeeze(dst_u_north)
                    dst_v_north = np.squeeze(dst_v_north)
                    dst_u_south = np.squeeze(dst_u_south)
                    dst_v_south = np.squeeze(dst_v_south)
                    dst_u_east = np.squeeze(dst_u_east)
                    dst_v_east = np.squeeze(dst_v_east)
                    dst_u_west = np.squeeze(dst_u_west)
                    dst_v_west = np.squeeze(dst_v_west)

                dst_ubar_north = np.zeros(dst_u_north.shape[1])
                dst_ubar_south = np.zeros(dst_u_south.shape[1])
                dst_ubar_east = np.zeros(dst_u_east.shape[1])
                dst_ubar_west = np.zeros(dst_u_west.shape[1])
                dst_vbar_north = np.zeros(dst_v_north.shape[1])
                dst_vbar_south = np.zeros(dst_v_south.shape[1])
                dst_vbar_east = np.zeros(dst_v_east.shape[1])
                dst_vbar_west = np.zeros(dst_v_west.shape[1])

#                print 'Shapes 3', dst_u_north.shape, dst_ubar_north.shape, z_u_north.shape, np.diff(z_u_north[:,1]).shape
                for i in range(dst_u_north.shape[1]):
                    dst_ubar_north[i] = (dst_u_north[:,i] * \
                        np.diff(z_u_north[:,i])).sum() / -z_u_north[0,i]
                    dst_ubar_south[i] = (dst_u_south[:,i] * \
                        np.diff(z_u_south[:,i])).sum() / -z_u_south[0,i]
                for i in range(dst_v_north.shape[1]):
                    dst_vbar_north[i] = (dst_v_north[:,i] * \
                        np.diff(z_v_north[:,i])).sum() / -z_v_north[0,i]
                    dst_vbar_south[i] = (dst_v_south[:,i] * \
                        np.diff(z_v_south[:,i])).sum() / -z_v_south[0,i]
                for j in range(dst_u_east.shape[1]):
                    dst_ubar_east[j] = (dst_u_east[:,j] * \
                        np.diff(z_u_east[:,j])).sum() / -z_u_east[0,j]
                    dst_ubar_west[j] = (dst_u_west[:,j] * \
                        np.diff(z_u_west[:,j])).sum() / -z_u_west[0,j]
                for j in range(dst_v_east.shape[1]):
                    dst_vbar_east[j] = (dst_v_east[:,j] * \
                        np.diff(z_v_east[:,j])).sum() / -z_v_east[0,j]
                    dst_vbar_west[j] = (dst_v_west[:,j] * \
                        np.diff(z_v_west[:,j])).sum() / -z_v_west[0,j]

                # spval
                idxu_north = np.where(dst_grd.hgrid.mask_u[-1,:] == 0)
                idxv_north = np.where(dst_grd.hgrid.mask_v[-1,:] == 0)
                idxu_south = np.where(dst_grd.hgrid.mask_u[0,:] == 0)
                idxv_south = np.where(dst_grd.hgrid.mask_v[0,:] == 0)
                idxu_east = np.where(dst_grd.hgrid.mask_u[:,-1] == 0)
                idxv_east = np.where(dst_grd.hgrid.mask_v[:,-1] == 0)
                idxu_west = np.where(dst_grd.hgrid.mask_u[:,0] == 0)
                idxv_west = np.where(dst_grd.hgrid.mask_v[:,0] == 0)

                dst_ubar_north[idxu_north[0]] = spval
                dst_vbar_north[idxv_north[0]] = spval
                dst_ubar_south[idxu_south[0]] = spval
                dst_vbar_south[idxv_south[0]] = spval
                dst_ubar_east[idxu_east[0]] = spval
                dst_vbar_east[idxv_east[0]] = spval
                dst_ubar_west[idxu_west[0]] = spval
                dst_vbar_west[idxv_west[0]] = spval

                nc.variables['ubar_north'][nctidx] = dst_ubar_north
                nc.variables['ubar_south'][nctidx] = dst_ubar_south
                nc.variables['ubar_east'][nctidx] = dst_ubar_east
                nc.variables['ubar_west'][nctidx] = dst_ubar_west

                nc.variables['vbar_north'][nctidx] = dst_vbar_north
                nc.variables['vbar_south'][nctidx] = dst_vbar_south
                nc.variables['vbar_east'][nctidx] = dst_vbar_east
                nc.variables['vbar_west'][nctidx] = dst_vbar_west

            nctidx = nctidx + 1
            print('ADDING to nctidx ', nctidx)
            nc.sync()
        # close files here? how?

    # close destination file
    nc.close()

    return
