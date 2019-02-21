import numpy as np
import os
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF
import matplotlib.pyplot as plt
import time
import datetime as dt
from matplotlib.dates import date2num, num2date

import pyroms
import pyroms_toolbox
import _remapping

class nctime(object):
    pass

def remap_uv(src_fileuv, src_grd, dst_grd, dst_fileu, dst_filev, dmax=0, cdepth=0, kk=0, dst_dir='./'):

    # CCS grid sub-sample
    xrange=src_grd.xrange; yrange=src_grd.yrange

    # get time
    nctime.long_name = 'time'
    nctime.units = 'days since 1900-01-01 00:00:00'

    # get dimensions
    Mp, Lp = dst_grd.hgrid.mask_rho.shape

    # create destination file
    print('\nCreating destination file', dst_fileu)
    if os.path.exists(dst_fileu) is True:
        os.remove(dst_fileu)
    pyroms_toolbox.nc_create_roms_file(dst_fileu, dst_grd, nctime)
    print('Creating destination file', dst_filev)
    if os.path.exists(dst_filev) is True:
        os.remove(dst_filev)
    pyroms_toolbox.nc_create_roms_file(dst_filev, dst_grd, nctime)

    # open destination file
    ncu = netCDF.Dataset(dst_fileu, 'a', format='NETCDF3_64BIT')
    ncv = netCDF.Dataset(dst_filev, 'a', format='NETCDF3_64BIT')

    #load var
    cdfuv = netCDF.Dataset(src_fileuv)
    src_varu = cdfuv.variables['u']
    src_varv = cdfuv.variables['v']

    tmp = cdfuv.variables['time'][:]
    if len(tmp) > 1:
        print('error : multiple frames in input file') ; exit()
    else:
        time = tmp[0]
    
    # we need to correct the time axis
    ref_soda = dt.datetime(1980,1,1,0,0)
    ref_roms = dt.datetime(1900,1,1,0,0)
    ndays = (ref_soda - ref_roms).days
    time = time + ndays 

    #get missing value
    spval = src_varu.missing_value

    # CCS grid sub-sample
    src_varu = src_varu[0,:, yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
    src_varv = src_varv[0,:, yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]

    # get weights file
    wts_file = 'remap_weights_' + src_grd.name + '_to_' + dst_grd.name + '_bilinear_uv_to_rho.nc'

    # build intermediate zgrid
    zlevel = -src_grd.z_t[::-1]
    nzlevel = len(zlevel)
    dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
    dst_grdz = pyroms.grid.ROMS_Grid(dst_grd.name+'_Z', dst_grd.hgrid, dst_zcoord)

    # create variable in destination file
    print('Creating variable u')
    ncu.createVariable('u', 'f8', ('ocean_time', 's_rho', 'eta_u', 'xi_u'), fill_value=spval)
    ncu.variables['u'].long_name = '3D u-momentum component'
    ncu.variables['u'].units = 'meter second-1'
    ncu.variables['u'].field = 'u-velocity, scalar, series'
    #ncu.variables['u_north']._FillValue = spval
    # create variable in destination file
    print('Creating variable ubar')
    ncu.createVariable('ubar', 'f8', ('ocean_time', 'eta_u', 'xi_u'), fill_value=spval)
    ncu.variables['ubar'].long_name = '2D u-momentum component'
    ncu.variables['ubar'].units = 'meter second-1'
    ncu.variables['ubar'].field = 'ubar-velocity,, scalar, series'
    #ncu.variables['ubar_north']._FillValue = spval

    print('Creating variable v')
    ncv.createVariable('v', 'f8', ('ocean_time', 's_rho', 'eta_v', 'xi_v'), fill_value=spval)
    ncv.variables['v'].long_name = '3D v-momentum component'
    ncv.variables['v'].units = 'meter second-1'
    ncv.variables['v'].field = 'v-velocity, scalar, series'
    #ncv.variables['v_north']._FillValue = spval
    print('Creating variable vbar')
    ncv.createVariable('vbar', 'f8', ('ocean_time', 'eta_v', 'xi_v'), fill_value=spval)
    ncv.variables['vbar'].long_name = '2D v-momentum component'
    ncv.variables['vbar'].units = 'meter second-1'
    ncv.variables['vbar'].field = 'vbar-velocity,, scalar, series'
    #ncv.variables['vbar_north']._FillValue = spval


    # remaping
    print('remapping and rotating u and v from', src_grd.name, \
                      'to', dst_grd.name)

    # flood the grid
    print('flood the grid')
    src_uz = pyroms_toolbox.BGrid_SODA.flood(src_varu, src_grd, Bpos='uv', \
                spval=spval, dmax=dmax, cdepth=cdepth, kk=kk)
    src_vz = pyroms_toolbox.BGrid_SODA.flood(src_varv, src_grd, Bpos='uv', \
                spval=spval, dmax=dmax, cdepth=cdepth, kk=kk)

    # horizontal interpolation using scrip weights
    print('horizontal interpolation using scrip weights')
    dst_uz = pyroms.remapping.remap(src_uz, wts_file, \
                                      spval=spval)
    dst_vz = pyroms.remapping.remap(src_vz, wts_file, \
                                      spval=spval)

    # vertical interpolation from standard z level to sigma
    print('vertical interpolation from standard z level to sigma')
    dst_u = pyroms.remapping.z2roms(dst_uz[::-1,:,:], dst_grdz, \
                        dst_grd, Cpos='rho', spval=spval, flood=False)
    dst_v = pyroms.remapping.z2roms(dst_vz[::-1,:,:], dst_grdz,  \
                        dst_grd, Cpos='rho', spval=spval, flood=False)


    # rotate u,v fields
    src_angle = np.zeros(dst_grd.hgrid.angle_rho.shape)
    dst_angle = dst_grd.hgrid.angle_rho
    angle = dst_angle - src_angle
    angle = np.tile(angle, (dst_grd.vgrid.N, 1, 1))
    U = dst_u + dst_v*1j
    eitheta = np.exp(-1j*angle[:,:,:])
    U = U * eitheta
    dst_u = np.real(U)
    dst_v = np.imag(U)


    # move back to u,v points
    dst_u = 0.5 * (dst_u[:,:,:-1] + dst_u[:,:,1:])
    dst_v = 0.5 * (dst_v[:,:-1,:] + dst_v[:,1:,:])

    # spval
    idxu = np.where(dst_grd.hgrid.mask_u == 0)
    idxv = np.where(dst_grd.hgrid.mask_v == 0)
    for n in range(dst_grd.vgrid.N):
        dst_u[n,idxu[0], idxu[1]] = spval
        dst_v[n,idxv[0], idxv[1]] = spval


    # compute depth average velocity ubar and vbar
    # get z at the right position
    z_u = 0.5 * (dst_grd.vgrid.z_w[0,:,:,:-1] + dst_grd.vgrid.z_w[0,:,:,1:])
    z_v = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,:] + dst_grd.vgrid.z_w[0,:,1:,:])

    dst_ubar = np.zeros((dst_u.shape[1], dst_u.shape[2]))
    dst_vbar = np.zeros((dst_v.shape[1], dst_v.shape[2]))

    for i in range(dst_ubar.shape[1]):
        for j in range(dst_ubar.shape[0]):
            dst_ubar[j,i] = (dst_u[:,j,i] * np.diff(z_u[:,j,i])).sum() / -z_u[0,j,i]

    for i in range(dst_vbar.shape[1]):
        for j in range(dst_vbar.shape[0]):
            dst_vbar[j,i] = (dst_v[:,j,i] * np.diff(z_v[:,j,i])).sum() / -z_v[0,j,i]

    # spval
    dst_ubar[idxu[0], idxu[1]] = spval
    dst_vbar[idxv[0], idxv[1]] = spval

    # write data in destination file
    print('write data in destination file\n')
    ncu.variables['ocean_time'][0] = time
    ncu.variables['u'][0] = dst_u
    ncu.variables['ubar'][0] = dst_ubar

    ncv.variables['ocean_time'][0] = time
    ncv.variables['v'][0] = dst_v
    ncv.variables['vbar'][0] = dst_vbar

    # close file
    ncu.close()
    ncv.close()
    cdfuv.close()
