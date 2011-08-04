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

def remapping(varname, srcfile, wts_files, srcgrd, dstgrd, \
              rotate_uv=False, trange=None, irange=None, jrange=None, \
              dstdir='./' ,zlevel=None, dmax=0, cdepth=0, kk=0, \
              uvar='u', vvar='v'):
    '''
    A remapping function to go from a ROMS grid to another ROMS grid.
    If the u/v variables need to be rotated, it must be called for each
    u/v pair (such as ubar/vbar, uice/vice).
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

    # build intermediaire zgrid
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
        raise ValueError, 'varname must be a str or a list of str'

    # if we're working on u and v, we'll compute ubar,vbar afterwards
    compute_ubar = False
    if varname.__contains__('u') == 1 and varname.__contains__('v') == 1:
        compute_ubar = True
        print 'ubar/vbar to be computed from u/v' 
        if varname.__contains__('ubar'):
            varname.remove('ubar')
            nvar = nvar-1
        if varname.__contains__('vbar'):
            varname.remove('vbar')
            nvar = nvar-1

    # if rotate_uv=True, check that u and v are in varname
    if rotate_uv is True:
        if varname.__contains__(uvar) == 0 or varname.__contains__(vvar) == 0:
            raise Warning, 'varname must include uvar and vvar in order to' \
                   + ' rotate the velocity field'
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
        raise ValueError, 'src_srcfile must be a str or a list of str'

    # get wts_file
    if type(wts_files).__name__ == 'str':
        wts_files = sorted(glob.glob(wts_files))
 
    # loop over the srcfile
    for nf in range(nfile):
        print 'Working with file', srcfile[nf], '...'

        # get time 
        ocean_time = pyroms.utility.get_nc_var('ocean_time', srcfile[nf])
        ntime = len(ocean_time[:])

        # trange argument
        if trange is None:
            trange = range(ntime)
        else:
            trange = range(trange[0], trange[1]+1)

        # create destination file
        dstfile = dstdir + os.path.basename(srcfile[nf])[:-3] + '_' + dstgrd.name + '.nc'
        if os.path.exists(dstfile) is False:
            print 'Creating destination file', dstfile
            pyroms_toolbox.nc_create_roms_file(dstfile, dstgrd, ocean_time)

        # open destination file
        nc = netCDF.Dataset(dstfile, 'a', format='NETCDF3_CLASSIC')

        nctidx = 0
        # loop over time
        for nt in trange:

            nc.variables['ocean_time'][nctidx] = ocean_time[nt]

            # loop over variable
            for nv in range(nvar):
                print ' '
                print 'remapping', varname[nv], 'from', srcgrd.name, \
                      'to', dstgrd.name
                print 'time =', ocean_time[nt]   

                # get source data
                src_var = pyroms.utility.get_nc_var(varname[nv], srcfile[nf])

                # determine variable dimension
                ndim = len(src_var.dimensions)-1

                # get spval
                try:
                    spval = src_var._FillValue
                except:
                    raise Warning, 'Did not find a _FillValue attribute.' 

                # irange
                if irange is None:
                    iirange = (0,src_var._shape()[-1])
                else:
                    iirange = irange

                # jrange
                if jrange is None:
                    jjrange = (0,src_var._shape()[-2])
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

                print 'Arakawa C-grid position is', Cpos

                # create variable in _destination file
                if nt == trange[0]:
                    print 'Creating variable', varname[nv]
                    nc.createVariable(varname[nv], 'f8', src_var.dimensions)
                    nc.variables[varname[nv]].long_name = src_var.long_name
                    try:
                        nc.variables[varname[nv]].units = src_var.units
                    except:
                        print varname[nv]+' has no units'
                    nc.variables[varname[nv]].time = src_var.time
                    nc.variables[varname[nv]].coordinates = \
                        src_var.coordinates
                    nc.variables[varname[nv]].field = src_var.field
                    nc.variables[varname[nv]]._FillValue = spval

                # get the right remap weights file
                for s in range(len(wts_files)):
                    if wts_files[s].__contains__(Cpos+'_to_'+Cpos+'.nc'):
                        wts_file = wts_files[s]
                        break
                    else:
                        if s == len(wts_files) - 1:
                            raise ValueError, 'Did not find the appropriate remap weights file'

                if ndim == 3:
                    # vertical interpolation from sigma to standard z level
                    print 'vertical interpolation from sigma to standard z level'
                    src_varz = pyroms.remapping.roms2z( \
                                 src_var[nt,:,jjrange[0]:jjrange[1],iirange[0]:iirange[1]], \
                                 srcgrd, srcgrdz, Cpos=Cpos, spval=spval, \
                                 irange=iirange, jrange=jjrange)

                    # flood the grid
                    print 'flood the grid'
                    src_varz = pyroms.remapping.flood(src_varz, srcgrdz, Cpos=Cpos, \
                                      irange=iirange, jrange=jjrange, spval=spval, \
                                      dmax=dmax, cdepth=cdepth, kk=kk)

                    tmp_src_varz = np.zeros((nzlevel, src_var[:].shape[-2], src_var[:].shape[-1]))
                    tmp_src_varz[:,jjrange[0]:jjrange[1],iirange[0]:iirange[1]] = src_varz
                else:
                    src_varz = src_var[nt,jjrange[0]:jjrange[1],iirange[0]:iirange[1]]
                    tmp_src_varz = np.zeros(src_var[nt,:].shape)
                    tmp_src_varz[jjrange[0]:jjrange[1],iirange[0]:iirange[1]] = src_varz

                print datetime.datetime.now()
                # horizontal interpolation using scrip weights
                print 'horizontal interpolation using scrip weights'
                dst_varz = pyroms.remapping.remap(tmp_src_varz, wts_file, \
                                                  spval=spval)

                print datetime.datetime.now()

                if ndim == 3:
                    # vertical interpolation from standard z level to sigma
                    print 'vertical interpolation from standard z level to sigma'
                    dst_var = pyroms.remapping.z2roms(dst_varz, dstgrdz, dstgrd, \
                                     Cpos=Cpos, spval=spval, flood=False)
                else:
                    dst_var = dst_varz

                # write data in destination file
                print 'write data in destination file'
                nc.variables[varname[nv]][nctidx] = dst_var
                if varname[nv] == 'u':
                    dst_u = dst_var
                if varname[nv] == 'v':
                    dst_v = dst_var


            # rotate the velocity field if requested
            if rotate_uv is True:
                print ' ' 
                print 'remapping and rotating u and v from', srcgrd.name, \
                      'to', dstgrd.name

                # get source data
                src_u = pyroms.utility.get_nc_var(uvar, srcfile[nf])
                src_v = pyroms.utility.get_nc_var(vvar, srcfile[nf])

                # get spval
                try:
                    spval = src_v._FillValue
                except:
                    raise Warning, 'Did not find a _FillValue attribute.' 

                # create variable in destination file
                print 'Creating variable '+uvar
                nc.createVariable(uvar, 'f8', src_u.dimensions)
                nc.variables[uvar].long_name = src_u.long_name
                nc.variables[uvar].units = src_u.units
                nc.variables[uvar].time = src_u.time
                nc.variables[uvar].coordinates = src_u.coordinates
                nc.variables[uvar].field = src_u.field
                nc.variables[uvar]._FillValue = spval
                print 'Creating variable '+vvar
                nc.createVariable(vvar, 'f8', src_v.dimensions)
                nc.variables[vvar].long_name = src_v.long_name
                nc.variables[vvar].units = src_v.units
                nc.variables[vvar].time = src_v.time
                nc.variables[uvar].coordinates = src_u.coordinates
                nc.variables[uvar].field = src_u.field
                nc.variables[vvar]._FillValue = spval

                # get the right remap weights file
                for s in range(len(wts_files)):
                    if wts_files[s].__contains__('u_to_rho.nc'):
                        wts_file_u = wts_files[s]
                    if wts_files[s].__contains__('v_to_rho.nc'):
                        wts_file_v = wts_files[s]

                # vertical interpolation from sigma to standard z level
                print 'vertical interpolation from sigma to standard z level'
                # irange
                if irange is None:
                    iirange = (0,src_u._shape()[-1])
                else:
                    iirange = irange

                # jrange
                if jrange is None:
                    jjrange = (0,src_u._shape()[-2])
                else:
                    jjrange = jrange

                ndim = len(src_v.dimensions)-1
                if ndim == 3:
                    src_uz = pyroms.remapping.roms2z( \
                            src_u[nt,:,jjrange[0]:jjrange[1],iirange[0]:iirange[1]], \
                            srcgrd, srcgrdz, Cpos='u', irange=iirange, jrange=jjrange)
                    # flood the grid
                    print 'flood the u grid'
                    src_uz = pyroms.remapping.flood(src_uz, srcgrdz, Cpos='u', \
                                  irange=iirange, jrange=jjrange, \
                                  spval=spval, dmax=dmax, cdepth=cdepth, kk=kk)
                    tmp_src_uz = np.zeros((nzlevel, src_u[:].shape[-2], src_u[:].shape[-1]))
                    tmp_src_uz[:,jjrange[0]:jjrange[1],iirange[0]:iirange[1]] = src_uz
                else:
                    src_uz = src_u[nt,jjrange[0]:jjrange[1],iirange[0]:iirange[1]]
                    tmp_src_uz = np.zeros(src_u[nt,:].shape)
                    tmp_src_uz[jjrange[0]:jjrange[1],iirange[0]:iirange[1]] = src_uz

                # irange
                if irange is None:
                    iirange = (0,src_v._shape()[-1])
                else:
                    iirange = irange

                # jrange
                if jrange is None:
                    jjrange = (0,src_v._shape()[-2])
                else:
                    jjrange = jrange

                if ndim == 3:
                    src_vz = pyroms.remapping.roms2z( \
                            src_v[nt,:,jjrange[0]:jjrange[1],iirange[0]:iirange[1]], \
                            srcgrd, srcgrdz, Cpos='v', irange=iirange, jrange=jjrange)

                    # flood the grid
                    print 'flood the v grid'
                    src_vz = pyroms.remapping.flood(src_vz, srcgrdz, Cpos='v', \
                                  irange=iirange, jrange=jjrange, \
                                  spval=spval, dmax=dmax, cdepth=cdepth, kk=kk)
                    tmp_src_vz = np.zeros((nzlevel, src_v[:].shape[-2], src_v[:].shape[-1]))
                    tmp_src_vz[:,jjrange[0]:jjrange[1],iirange[0]:iirange[1]] = src_vz
                else:
                    src_vz = src_v[nt,jjrange[0]:jjrange[1],iirange[0]:iirange[1]]
                    tmp_src_vz = np.zeros(src_v[nt,:].shape)
                    tmp_src_vz[jjrange[0]:jjrange[1],iirange[0]:iirange[1]] = src_vz

                # horizontal interpolation using scrip weights
                print 'horizontal interpolation using scrip weights'
                dst_uz = pyroms.remapping.remap(tmp_src_uz, wts_file_u, \
                                                  spval=spval)
                dst_vz = pyroms.remapping.remap(tmp_src_vz, wts_file_v, \
                                                  spval=spval)

                if ndim == 3:
                    # vertical interpolation from standard z level to sigma
                    print 'vertical interpolation from standard z level to sigma'
                    dst_u = pyroms.remapping.z2roms(dst_uz, dstgrdz, dstgrd, \
                                 Cpos='rho', spval=spval, flood=False)
                    dst_v = pyroms.remapping.z2roms(dst_vz, dstgrdz, dstgrd, \
                                 Cpos='rho', spval=spval, flood=False)
                else:
                    dst_u = dst_uz
                    dst_v = dst_vz

                # rotate u,v fields
                for s in range(len(wts_files)):
                    if wts_files[s].__contains__('rho_to_rho.nc'):
                        wts_file = wts_files[s]
                src_angle = pyroms.remapping.remap(srcgrd.hgrid.angle_rho, wts_file)
                dst_angle = dstgrd.hgrid.angle_rho
                angle = dst_angle - src_angle
                angle = np.tile(angle, (dstgrd.vgrid.N, 1, 1))

                U = dst_u + dst_v*1j
                eitheta = np.exp(-1j*angle)
                U = U * eitheta

                dst_u = np.real(U)
                dst_v = np.imag(U)

                # move back to u,v points
                dst_u = 0.5 * (dst_u[:,:,:-1] + dst_u[:,:,1:])
                dst_v = 0.5 * (dst_v[:,:-1,:] + dst_v[:,1:,:])

                # spval
                idxu = np.where(dstgrd.hgrid.mask_u == 0)
                idxv = np.where(dstgrd.hgrid.mask_v == 0)
                for n in range(dstgrd.vgrid.N):
                    dst_u[n,idxu[0], idxu[1]] = spval
                    dst_v[n,idxv[0], idxv[1]] = spval

                # write data in destination file
                print 'write data in destination file'
                nc.variables[uvar][nctidx] = dst_u
                nc.variables[vvar][nctidx] = dst_v

            if compute_ubar:
                print 'Creating variable ubar'
                nc.createVariable('ubar', 'f8', \
                     ('ocean_time', 'eta_u', 'xi_u'), fill_value=spval)
                nc.variables['ubar'].long_name = '2D u-momentum component'
                nc.variables['ubar'].units = 'meter second-1'
                nc.variables['ubar'].time = 'ocean_time'
                nc.variables['ubar'].coordinates = 'xi_u eta_u ocean_time'
                nc.variables['ubar'].field = 'ubar-velocity,, scalar, series'
                nc.variables['ubar']._FillValue = spval
                print 'Creating variable vbar'
                nc.createVariable('vbar', 'f8', \
                     ('ocean_time', 'eta_v', 'xi_v'), fill_value=spval)
                nc.variables['vbar'].long_name = '2D v-momentum component'
                nc.variables['vbar'].units = 'meter second-1'
                nc.variables['vbar'].time = 'ocean_time'
                nc.variables['vbar'].coordinates = 'xi_v eta_v ocean_time'
                nc.variables['vbar'].field = 'vbar-velocity,, scalar, series'
                nc.variables['vbar']._FillValue = spval

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

                nc.variables['ubar'][0] = dst_ubar
                nc.variables['vbar'][0] = dst_vbar

        nctidx = nctidx + 1
 
    # close destination file
    nc.close()

    return
