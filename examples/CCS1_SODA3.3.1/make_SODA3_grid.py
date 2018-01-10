#!/usr/bin/env python

import netCDF4 as nc
import numpy as np
import matplotlib.pylab as plt

# use any output file for lon/lat
filein='/Volumes/P1/Data/SODA/SODA_3.3.1/soda3.3.1_5dy_ocean_reg_1980_01_08.nc'

fid = nc.Dataset(filein,'r')
lon_t = fid.variables['xt_ocean'][:]
lat_t = fid.variables['yt_ocean'][:]
lon_uv = fid.variables['xu_ocean'][:]
lat_uv = fid.variables['yu_ocean'][:]
st_ocean = fid.variables['st_ocean'][:]
temp = fid.variables['temp'][:].squeeze()
u = fid.variables['u'][:].squeeze()
fid.close()

#--- dimensions
nz = st_ocean.shape[0]
ny = lat_t.shape[0]
nx = lon_t.shape[0]
spval=-1.0e+20

#--- create 2d lon/lat arrays
lon_t_2d, lat_t_2d = np.meshgrid(lon_t,lat_t)
lon_uv_2d, lat_uv_2d = np.meshgrid(lon_uv,lat_uv)

#--- coriolis
deg2rad = 2 * np.pi / 360
omega = 2 * np.pi / 86400
coriolis_param = omega * np.sin(deg2rad * lat_t_2d)

#--- create depth arrays
sw_ocean = np.empty((nz))
sw_ocean[:-1] = 0.5 * (st_ocean[:-1] + st_ocean[1:])
sw_ocean[-1] = 5500.

st_edges_ocean = np.empty((nz+1))
st_edges_ocean[1:] = sw_ocean[:]
st_edges_ocean[0] = 0.

sw_edges_ocean = np.empty((nz+1))
sw_edges_ocean[:-1] = st_ocean[:]
sw_edges_ocean[-1] = 5500.

#--- compute kmt,...

mask_t = np.ones((nz,ny,nx))
mask_t[np.where(temp.mask == True)] = 0

mask_u = np.ones((nz,ny,nx))
mask_u[np.where(u.mask == True)] = 0

kmt = np.empty((ny,nx))
ht = np.empty((ny,nx))
for ky in np.arange(ny):
    for kx in np.arange(nx):
        indlist = np.where(mask_t[:,ky,kx] == 1)[0]
        if len(indlist) == 0:
            kmt[ky,kx] = 0
            ht[ky,kx] = 0
        else:
            kmt[ky,kx] = indlist.max() + 1
            ht[ky,kx] = st_edges_ocean[indlist.max() + 1]

kmt[np.where(kmt == 0)] = spval
ht[np.where(ht == 0)] = spval

kmu = np.empty((ny,nx))
for ky in np.arange(ny):
    for kx in np.arange(nx):
        indlist = np.where(mask_u[:,ky,kx] == 1)[0]
        if len(indlist) == 0:
            kmu[ky,kx] = 0
        else:
            kmu[ky,kx] = indlist.max() + 1

kmu[np.where(kmu == 0)] = spval

#--- write the output file
fidout = nc.Dataset('SODA3_0.5deg_grid.nc','w',format='NETCDF3_CLASSIC')
fidout.createDimension('time', None)
fidout.createDimension('xu_ocean',nx)
fidout.createDimension('xt_ocean',nx)
fidout.createDimension('yu_ocean',ny)
fidout.createDimension('yt_ocean',ny)
fidout.createDimension('st_ocean',nz)
fidout.createDimension('sw_ocean',nz)
fidout.createDimension('st_edges_ocean',nz+1)
fidout.createDimension('sw_edges_ocean',nz+1)

nc_lon_t = fidout.createVariable('geolon_t','f8',('yt_ocean','xt_ocean',))
nc_lat_t = fidout.createVariable('geolat_t','f8',('yt_ocean','xt_ocean',))
nc_lon_uv = fidout.createVariable('geolon_c','f8',('yu_ocean','xu_ocean',))
nc_lat_uv = fidout.createVariable('geolat_c','f8',('yu_ocean','xu_ocean',))
nc_st_ocean = fidout.createVariable('st_ocean','f8',('st_ocean',))
nc_sw_ocean = fidout.createVariable('sw_ocean','f8',('sw_ocean',))
nc_st_edges_ocean = fidout.createVariable('st_edges_ocean','f8',('st_edges_ocean',))
nc_sw_edges_ocean = fidout.createVariable('sw_edges_ocean','f8',('sw_edges_ocean',))
nc_coriolis_param = fidout.createVariable('coriolis_param','f8',('yt_ocean','xt_ocean',))
nc_kmt = fidout.createVariable('kmt','f8',('yt_ocean','xt_ocean',),fill_value=spval)
nc_ht = fidout.createVariable('ht','f8',('yt_ocean','xt_ocean',),fill_value=spval)
nc_kmu = fidout.createVariable('kmu','f8',('yu_ocean','xu_ocean',),fill_value=spval)

nc_lon_t[:]  = lon_t_2d
nc_lat_t[:]  = lat_t_2d
nc_lon_uv[:] = lon_uv_2d
nc_lat_uv[:] = lat_uv_2d

nc_st_ocean[:]  = st_ocean
nc_sw_ocean[:]  = sw_ocean
nc_st_edges_ocean[:]  = st_edges_ocean
nc_sw_edges_ocean[:]  = sw_edges_ocean

nc_coriolis_param[:]  = coriolis_param
nc_kmt[:]             = kmt
nc_ht[:]             = ht
nc_kmu[:]             = kmu

fidout.close()
