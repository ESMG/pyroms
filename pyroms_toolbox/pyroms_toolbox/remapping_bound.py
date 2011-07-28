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

def remapping_bound(varname, srcfile, wts_files, srcgrd, dstgrd, \
              rotate_uv=False, trange=None, irange=None, jrange=None, \
              dstdir='./' ,zlevel=None, dmax=0, cdepth=0, kk=0, \
              uvar='u', vvar='v'):
    '''
    A remapping function to extract boundary conditions from one ROMS grid
    to another. It will optionally rotating u and v variables, but needs
    to be called separately for each u/v pair (such as ubar/vbar, uice/vice).
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
        dstfile = dstdir + os.path.basename(srcfile[nf])[:-3] + '_' \
               + dstgrd.name + '_bdry.nc'
        if os.path.exists(dstfile) is False:
            print 'Creating destination file', dstfile
            pyroms_toolbox.nc_create_roms_file(dstfile, dstgrd, \
                ocean_time, lgrid=False)

        # open destination file
        nc = netCDF.Dataset(dstfile, 'a', format='NETCDF3_CLASSIC')

        sides = ['_west','_east','_north','_south']
	long = {'_west':'Western', '_east':'Eastern', \
	        '_north':'Northern', '_south':'Southern'}
	dimexcl = {'_west':'xi', '_east':'xi', \
	           '_north':'eta', '_south':'eta'}
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
		    for sid in sides:
		       varn = varname[nv]+str[sid]
                       print 'Creating variable', varn
		       dims = src_var.dimensions
		       for dim in dims:
		           if dim.match(dimexcl[sid]):
			       dims.remove(dim)
                       nc.createVariable(varn, 'f8', dims)
                       nc.variables[varn].long_name = varname[nv] + \
		            ' ' + long[sid] + ' boundary condition'
                       try:
                           nc.variables[varn].units = src_var.units
                       except:
                           print varn+' has no units'
                       nc.variables[varn].time = src_var.time
                       nc.variables[varn].coordinates = \
                           str(dims.reverse())
                       nc.variables[varn].field = src_var.field
                       nc.variables[varn]._FillValue = spval

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

                if ndim == 3:
                    dst_var_north = pyroms.remapping.z2roms(dst_varz[::-1, \
                          Mp-1:Mp,0:Lp], dst_grdz, dst_grd, Cpos=Cpos, \
                          spval=spval, flood=False, irange=(0,Lp), \
                          jrange=(Mp-1,Mp))
                    dst_var_south = pyroms.remapping.z2roms(dst_varz[::-1, \
                          0:1, :], dst_grdz, dst_grd, Cpos=Cpos, \
                          spval=spval, flood=False, irange=(0,Lp), \
                          jrange=(0,1))
                    dst_var_east = pyroms.remapping.z2roms(dst_varz[::-1, \
                          :, Lp-1:Lp], dst_grdz, dst_grd, Cpos=Cpos, \
                          spval=spval, flood=False, irange=(Lp-1,Lp), \
                          jrange=(0,Mp))
                    dst_var_west = pyroms.remapping.z2roms(dst_varz[::-1, \
                          :, 0:1], dst_grdz, dst_grd, Cpos=Cpos, \
                          spval=spval, flood=False, irange=(0,1), \
                          jrange=(0,Mp))
                else:
                    dst_var_north = dst_varz[-1, :]
                    dst_var_south = dst_varz[0, :]
                    dst_var_east = dst_varz[:, -1]
                    dst_var_west = dst_varz[:, 0]

                print datetime.datetime.now()

                # write data in destination file
                print 'write data in destination file'
		sid = '_west'
		varn = varname[nv]+str[sid]
                nc.variables[varn][nctidx] = dst_var_west

		sid = '_east'
		varn = varname[nv]+str[sid]
                nc.variables[varn][nctidx] = dst_var_east

		sid = '_north'
		varn = varname[nv]+str[sid]
                nc.variables[varn][nctidx] = dst_var_north

		sid = '_south'
		varn = varname[nv]+str[sid]
                nc.variables[varn][nctidx] = dst_var_south

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

                # write data in destination file
                print 'write data in destination file'
                nc.variables[uvar][nctidx] = dst_u
                nc.variables[vvar][nctidx] = dst_v

        nctidx = nctidx + 1
 
    # close destination file
    nc.close()

    return
