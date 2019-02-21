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

def remap_bdry_uv(src_fileuv, src_grd, dst_grd, dst_fileu, dst_filev, dmax=0, cdepth=0, kk=0, dst_dir='./'):

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
    pyroms_toolbox.nc_create_roms_bdry_file(dst_filev, dst_grd, nctime)

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
    print('Creating variable u_north')
    ncu.createVariable('u_north', 'f8', ('ocean_time', 's_rho', 'xi_u'), fill_value=spval)
    ncu.variables['u_north'].long_name = '3D u-momentum north boundary condition'
    ncu.variables['u_north'].units = 'meter second-1'
    ncu.variables['u_north'].field = 'u_north, scalar, series'
    #ncu.variables['u_north']._FillValue = spval
    print('Creating variable u_south')
    ncu.createVariable('u_south', 'f8', ('ocean_time', 's_rho', 'xi_u'), fill_value=spval)
    ncu.variables['u_south'].long_name = '3D u-momentum south boundary condition'
    ncu.variables['u_south'].units = 'meter second-1'
    ncu.variables['u_south'].field = 'u_south, scalar, series'
    #ncu.variables['u_south']._FillValue = spval
    print('Creating variable u_west')
    ncu.createVariable('u_west', 'f8', ('ocean_time', 's_rho', 'eta_u'), fill_value=spval)
    ncu.variables['u_west'].long_name = '3D u-momentum west boundary condition'
    ncu.variables['u_west'].units = 'meter second-1'
    ncu.variables['u_west'].field = 'u_east, scalar, series'
    #ncu.variables['u_west']._FillValue = spval

    # create variable in destination file
    print('Creating variable ubar_north')
    ncu.createVariable('ubar_north', 'f8', ('ocean_time', 'xi_u'), fill_value=spval)
    ncu.variables['ubar_north'].long_name = '2D u-momentum north boundary condition'
    ncu.variables['ubar_north'].units = 'meter second-1'
    ncu.variables['ubar_north'].field = 'ubar_north, scalar, series'
    #ncu.variables['ubar_north']._FillValue = spval
    print('Creating variable ubar_south')
    ncu.createVariable('ubar_south', 'f8', ('ocean_time', 'xi_u'), fill_value=spval)
    ncu.variables['ubar_south'].long_name = '2D u-momentum south boundary condition'
    ncu.variables['ubar_south'].units = 'meter second-1'
    ncu.variables['ubar_south'].field = 'ubar_south, scalar, series'
    #ncu.variables['ubar_south']._FillValue = spval
    print('Creating variable ubar_west')
    ncu.createVariable('ubar_west', 'f8', ('ocean_time', 'eta_u'), fill_value=spval)
    ncu.variables['ubar_west'].long_name = '2D u-momentum west boundary condition'
    ncu.variables['ubar_west'].units = 'meter second-1'
    ncu.variables['ubar_west'].field = 'ubar_east, scalar, series'
    #ncu.variables['ubar_west']._FillValue = spval

    print('Creating variable v_north')
    ncv.createVariable('v_north', 'f8', ('ocean_time', 's_rho', 'xi_v'), fill_value=spval)
    ncv.variables['v_north'].long_name = '3D v-momentum north boundary condition'
    ncv.variables['v_north'].units = 'meter second-1'
    ncv.variables['v_north'].field = 'v_north, scalar, series'
    #ncv.variables['v_north']._FillValue = spval
    print('Creating variable v_south')
    ncv.createVariable('v_south', 'f8', ('ocean_time', 's_rho', 'xi_v'), fill_value=spval)
    ncv.variables['v_south'].long_name = '3D v-momentum south boundary condition'
    ncv.variables['v_south'].units = 'meter second-1'
    ncv.variables['v_south'].field = 'v_south, scalar, series'
    #ncv.variables['v_south']._FillValue = spval
    print('Creating variable v_west')
    ncv.createVariable('v_west', 'f8', ('ocean_time', 's_rho', 'eta_v'), fill_value=spval)
    ncv.variables['v_west'].long_name = '3D v-momentum west boundary condition'
    ncv.variables['v_west'].units = 'meter second-1'
    ncv.variables['v_west'].field = 'v_east, scalar, series'
    #ncv.variables['v_west']._FillValue = spval

    print('Creating variable vbar_north')
    ncv.createVariable('vbar_north', 'f8', ('ocean_time', 'xi_v'), fill_value=spval)
    ncv.variables['vbar_north'].long_name = '2D v-momentum north boundary condition'
    ncv.variables['vbar_north'].units = 'meter second-1'
    ncv.variables['vbar_north'].field = 'vbar_north, scalar, series'
    #ncv.variables['vbar_north']._FillValue = spval
    print('Creating variable vbar_south')
    ncv.createVariable('vbar_south', 'f8', ('ocean_time', 'xi_v'), fill_value=spval)
    ncv.variables['vbar_south'].long_name = '2D v-momentum south boundary condition'
    ncv.variables['vbar_south'].units = 'meter second-1'
    ncv.variables['vbar_south'].field = 'vbar_south, scalar, series'
    #ncv.variables['vbar_south']._FillValue = spval
    print('Creating variable vbar_west')
    ncv.createVariable('vbar_west', 'f8', ('ocean_time', 'eta_v'), fill_value=spval)
    ncv.variables['vbar_west'].long_name = '2D v-momentum west boundary condition'
    ncv.variables['vbar_west'].units = 'meter second-1'
    ncv.variables['vbar_west'].field = 'vbar_east, scalar, series'
    #ncv.variables['vbar_west']._FillValue = spval



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
    dst_u_north = pyroms.remapping.z2roms(dst_uz[::-1, Mp-2:Mp, 0:Lp], \
                      dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                      flood=False, irange=(0,Lp), jrange=(Mp-2,Mp))
    dst_u_south = pyroms.remapping.z2roms(dst_uz[::-1, 0:2, 0:Lp], \
                      dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                      flood=False, irange=(0,Lp), jrange=(0,2))
    dst_u_west = pyroms.remapping.z2roms(dst_uz[::-1, 0:Mp, 0:2], \
                      dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                      flood=False, irange=(0,2), jrange=(0,Mp))

    dst_v_north = pyroms.remapping.z2roms(dst_vz[::-1, Mp-2:Mp, 0:Lp], \
                      dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                      flood=False, irange=(0,Lp), jrange=(Mp-2,Mp))
    dst_v_south = pyroms.remapping.z2roms(dst_vz[::-1, 0:2, 0:Lp], \
                      dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                      flood=False, irange=(0,Lp), jrange=(0,2))
    dst_v_west = pyroms.remapping.z2roms(dst_vz[::-1, 0:Mp, 0:2], \
                      dst_grdz, dst_grd, Cpos='rho', spval=spval, \
                      flood=False, irange=(0,2), jrange=(0,Mp))

    # rotate u,v fields
    src_angle = np.zeros(dst_grd.hgrid.angle_rho.shape)
    dst_angle = dst_grd.hgrid.angle_rho
    angle = dst_angle - src_angle
    angle = np.tile(angle, (dst_grd.vgrid.N, 1, 1))

    U_north = dst_u_north + dst_v_north*1j
    eitheta_north = np.exp(-1j*angle[:,Mp-2:Mp, 0:Lp])
    U_north = U_north * eitheta_north
    dst_u_north = np.real(U_north)
    dst_v_north = np.imag(U_north)

    U_south = dst_u_south + dst_v_south*1j
    eitheta_south = np.exp(-1j*angle[:,0:2, 0:Lp])
    U_south = U_south * eitheta_south
    dst_u_south = np.real(U_south)
    dst_v_south = np.imag(U_south)

    U_west = dst_u_west + dst_v_west*1j
    eitheta_west = np.exp(-1j*angle[:,0:Mp, 0:2])
    U_west = U_west * eitheta_west
    dst_u_west = np.real(U_west)
    dst_v_west = np.imag(U_west)

    # move back to u,v points
    dst_u_north = 0.5 * np.squeeze(dst_u_north[:,-1,:-1] + dst_u_north[:,-1,1:])
    dst_v_north = 0.5 * np.squeeze(dst_v_north[:,:-1,:] + dst_v_north[:,1:,:])
    dst_u_south = 0.5 * np.squeeze(dst_u_south[:,0,:-1] + dst_u_south[:,0,1:])
    dst_v_south = 0.5 * np.squeeze(dst_v_south[:,:-1,:] + dst_v_south[:,1:,:])
    dst_u_west = 0.5 * np.squeeze(dst_u_west[:,:,:-1] + dst_u_west[:,:,1:])
    dst_v_west = 0.5 * np.squeeze(dst_v_west[:,:-1,0] + dst_v_west[:,1:,0])

    # spval
    idxu_north = np.where(dst_grd.hgrid.mask_u[-1,:] == 0)
    idxv_north = np.where(dst_grd.hgrid.mask_v[-1,:] == 0)
    idxu_south = np.where(dst_grd.hgrid.mask_u[0,:] == 0)
    idxv_south = np.where(dst_grd.hgrid.mask_v[0,:] == 0)
    idxu_west = np.where(dst_grd.hgrid.mask_u[:,0] == 0)
    idxv_west = np.where(dst_grd.hgrid.mask_v[:,0] == 0)
    for n in range(dst_grd.vgrid.N):
        dst_u_north[n, idxu_north[0]] = spval
        dst_v_north[n, idxv_north[0]] = spval
        dst_u_south[n, idxu_south[0]] = spval
        dst_v_south[n, idxv_south[0]] = spval
        dst_u_west[n, idxu_west[0]] = spval
        dst_v_west[n, idxv_west[0]] = spval

    # compute depth average velocity ubar and vbar
    # get z at the right position
    z_u_north = 0.5 * (dst_grd.vgrid.z_w[0,:,-1,:-1] + dst_grd.vgrid.z_w[0,:,-1,1:])
    z_v_north = 0.5 * (dst_grd.vgrid.z_w[0,:,-1,:] + dst_grd.vgrid.z_w[0,:,-2,:])
    z_u_south = 0.5 * (dst_grd.vgrid.z_w[0,:,0,:-1] + dst_grd.vgrid.z_w[0,:,0,1:])
    z_v_south = 0.5 * (dst_grd.vgrid.z_w[0,:,0,:] + dst_grd.vgrid.z_w[0,:,1,:])
    z_u_west = 0.5 * (dst_grd.vgrid.z_w[0,:,:,0] + dst_grd.vgrid.z_w[0,:,:,1])
    z_v_west = 0.5 * (dst_grd.vgrid.z_w[0,:,:-1,0] + dst_grd.vgrid.z_w[0,:,1:,0])

    dst_ubar_north = np.zeros(dst_u_north.shape[1])
    dst_ubar_south = np.zeros(dst_u_south.shape[1])
    dst_ubar_west = np.zeros(dst_u_west.shape[1])
    dst_vbar_north = np.zeros(dst_v_north.shape[1])
    dst_vbar_south = np.zeros(dst_v_south.shape[1])
    dst_vbar_west = np.zeros(dst_v_west.shape[1])

    for i in range(dst_u_north.shape[1]):
        dst_ubar_north[i] = (dst_u_north[:,i] * np.diff(z_u_north[:,i])).sum() / -z_u_north[0,i]
        dst_ubar_south[i] = (dst_u_south[:,i] * np.diff(z_u_south[:,i])).sum() / -z_u_south[0,i]
    for i in range(dst_v_north.shape[1]):
        dst_vbar_north[i] = (dst_v_north[:,i] * np.diff(z_v_north[:,i])).sum() / -z_v_north[0,i]
        dst_vbar_south[i] = (dst_v_south[:,i] * np.diff(z_v_south[:,i])).sum() / -z_v_south[0,i]
    for j in range(dst_u_west.shape[1]):
        dst_ubar_west[j] = (dst_u_west[:,j] * np.diff(z_u_west[:,j])).sum() / -z_u_west[0,j]
    for j in range(dst_v_west.shape[1]):
        dst_vbar_west[j] = (dst_v_west[:,j] * np.diff(z_v_west[:,j])).sum() / -z_v_west[0,j]

    #mask
    dst_ubar_north = np.ma.masked_where(dst_grd.hgrid.mask_u[-1,:] == 0, dst_ubar_north)
    dst_ubar_south = np.ma.masked_where(dst_grd.hgrid.mask_u[0,:] == 0, dst_ubar_south)
    dst_ubar_west = np.ma.masked_where(dst_grd.hgrid.mask_u[:,0] == 0, dst_ubar_west)
    dst_vbar_north = np.ma.masked_where(dst_grd.hgrid.mask_v[-1,:] == 0, dst_vbar_north)
    dst_vbar_south = np.ma.masked_where(dst_grd.hgrid.mask_v[0,:] == 0, dst_vbar_south)
    dst_vbar_west = np.ma.masked_where(dst_grd.hgrid.mask_v[:,0] == 0, dst_vbar_west)

    # write data in destination file
    print('write data in destination file')
    ncu.variables['ocean_time'][0] = time
    ncu.variables['u_north'][0] = dst_u_north
    ncu.variables['u_south'][0] = dst_u_south
    ncu.variables['u_west'][0] = dst_u_west
    ncu.variables['ubar_north'][0] = dst_ubar_north
    ncu.variables['ubar_south'][0] = dst_ubar_south
    ncu.variables['ubar_west'][0] = dst_ubar_west

    ncv.variables['ocean_time'][0] = time
    ncv.variables['v_north'][0] = dst_v_north
    ncv.variables['v_south'][0] = dst_v_south
    ncv.variables['v_west'][0] = dst_v_west
    ncv.variables['vbar_north'][0] = dst_vbar_north
    ncv.variables['vbar_south'][0] = dst_vbar_south
    ncv.variables['vbar_west'][0] = dst_vbar_west

    # close file
    ncu.close()
    ncv.close()
    cdfuv.close()
