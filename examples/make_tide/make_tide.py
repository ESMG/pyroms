'''
this is an example to generate a ROMS tide file for TPXO8-atlas tidal product. Till this moment I couldn't found
a good python script to do so, so this script is created for this purpose. TPXO8-atlas is a high resolution global
model of ocean tides. it is available here:

http://volkov.oce.orst.edu/tides/tpxo8_atlas.html

download the 8 major constitutes with 1/30 deg resotuion. Other constitutes are also supported, but haven't been
tested. Only netCDF format is supported in this version, so make sure you download the right files.

The script 'make_remap_weights_file.py' is used to generate weight file between two grids, in this example the
CHUKCHI2 and TPXO8 grid files. TPXO8 grid file is kept in the same directory as data files:

/archive/u1/uaf/kate/tides/tpxo8/tpxo8_atlas_30c_v1.nc

Run this script before 'make_tide.py' to generate the weight files (remap_*.nc).

To use the remap function a new set of grid is created (CGrid_TPXO8). I followed the structure of GLORY CGrid
from the PYROMS packaage itself, so that it is consistent with PYROMS.

---------------------------------------------------------------------------------------------------------------------
tidal parameters are converted to elliptical parameters using ap2ep, origionally written in MATLAB and later
translated to python package PyFVCOM. See

CGrid_TPXO8/tidal_ellipse.py for more details.

-Chuning Wang, Dec 18 2016

'''

import netCDF4 as nc
import numpy as np
from scipy.interpolate import griddata

import pyroms
import pyroms_toolbox

import CGrid_TPXO8

# tidal constituents
consts1 =['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2']
consts2 = ['MF']
consts = consts1+consts2
# consts =['M2']
consts_num = len(consts)

# read ROMS grid
dstgrd = pyroms.grid.get_ROMS_grid('CHUKCHI2')
lat = dstgrd.hgrid.lat_rho
lon = dstgrd.hgrid.lon_rho
msk = dstgrd.hgrid.mask_rho
eta, xi = lat.shape

# read weights file
wts_file_t = 'remap_weights_TPXO8_to_CHUKCHI2_bilinear_t_to_rho.nc'
wts_file_u = 'remap_weights_TPXO8_to_CHUKCHI2_bilinear_u_to_rho.nc'
wts_file_v = 'remap_weights_TPXO8_to_CHUKCHI2_bilinear_v_to_rho.nc'

wts_file_t_lr = 'remap_weights_TPXO8atlas6_to_CHUKCHI2_bilinear_t_to_rho.nc'
wts_file_u_lr = 'remap_weights_TPXO8atlas6_to_CHUKCHI2_bilinear_u_to_rho.nc'
wts_file_v_lr = 'remap_weights_TPXO8atlas6_to_CHUKCHI2_bilinear_v_to_rho.nc'

# read TPXO8 files
pth_tpxo = '/archive/u1/uaf/kate/tides/tpxo8/'
srcgrd = CGrid_TPXO8.get_nc_CGrid_TPXO8(pth_tpxo+'grid_tpxo8atlas_30_v1.nc', \
      xrange=(4400, 5100), yrange=(5200, 6500))
srcgrd_lr = CGrid_TPXO8.get_nc_CGrid_TPXO8(pth_tpxo+'grid_tpxo8_atlas6.nc', name='TPXO8atlas6', \
      xrange=(920, 1020), yrange=(1030, 1300))

missing_value = srcgrd.missing_value

# grid sub-sample
xrange = srcgrd.xrange
yrange = srcgrd.yrange

xrange_lr = srcgrd_lr.xrange
yrange_lr = srcgrd_lr.yrange

# initiate variables
hamp = np.zeros((consts_num, eta, xi))*np.nan
hpha = np.zeros((consts_num, eta, xi))*np.nan
cmax = np.zeros((consts_num, eta, xi))*np.nan
cmin = np.zeros((consts_num, eta, xi))*np.nan
cang = np.zeros((consts_num, eta, xi))*np.nan
cpha = np.zeros((consts_num, eta, xi))*np.nan

k = -1
for cst in consts:
    k = k+1

    if cst in consts1:
        fid = nc.Dataset(pth_tpxo+'hf.'+cst.lower()+'_tpxo8_atlas_30c_v1.nc', 'r')
        hRe = fid.variables['hRe'][yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        hIm = fid.variables['hIm'][yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
    elif cst in consts2:
        fid = nc.Dataset(pth_tpxo+'hf.'+cst.lower()+'_tpxo8_atlas_6.nc', 'r')
        hRe = fid.variables['hRe'][yrange_lr[0]:yrange_lr[1]+1, xrange_lr[0]:xrange_lr[1]+1]
        hIm = fid.variables['hIm'][yrange_lr[0]:yrange_lr[1]+1, xrange_lr[0]:xrange_lr[1]+1]
    fid.close()

    if cst in consts1:
        fid = nc.Dataset(pth_tpxo+'uv.'+cst.lower()+'_tpxo8_atlas_30c_v1.nc', 'r')
        uRe = fid.variables['uRe'][yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        uIm = fid.variables['uIm'][yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        vRe = fid.variables['vRe'][yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
        vIm = fid.variables['vIm'][yrange[0]:yrange[1]+1, xrange[0]:xrange[1]+1]
    elif cst in consts2:
        fid = nc.Dataset(pth_tpxo+'uv.'+cst.lower()+'_tpxo8_atlas_6.nc', 'r')
        uRe = fid.variables['uRe'][yrange_lr[0]:yrange_lr[1]+1, xrange_lr[0]:xrange_lr[1]+1]
        uIm = fid.variables['uIm'][yrange_lr[0]:yrange_lr[1]+1, xrange_lr[0]:xrange_lr[1]+1]
        vRe = fid.variables['vRe'][yrange_lr[0]:yrange_lr[1]+1, xrange_lr[0]:xrange_lr[1]+1]
        vIm = fid.variables['vIm'][yrange_lr[0]:yrange_lr[1]+1, xrange_lr[0]:xrange_lr[1]+1]
    fid.close()

    hRe = hRe*1.0
    hIm = hIm*1.0
    uRe = uRe*1.0
    uIm = uIm*1.0
    vRe = vRe*1.0
    vIm = vIm*1.0

    # unit convertion
    hRe = hRe/1000.  # mm to m
    hIm = hIm/1000.  # mm to m
    uRe = uRe/10000.  # cm2s-1 to m2s-1
    uIm = uIm/10000.  # cm2s-1 to m2s-1
    vRe = vRe/10000.  # cm2s-1 to m2s-1
    vIm = vIm/10000.  # cm2s-1 to m2s-1

    # convert vertical transport (uRe, vRe) to velocity
    if cst in consts1:
        zu = srcgrd.z_u
        zv = srcgrd.z_v
    elif cst in consts2:
        zu = srcgrd_lr.z_u
        zv = srcgrd_lr.z_v

    uRe = uRe/zu  # ms-1
    uIm = uIm/zu  # ms-1
    vRe = vRe/zv  # ms-1
    vIm = vIm/zv  # ms-1

    if cst in consts1:
        # mask invalid data
        hRe = np.ma.masked_where(~srcgrd.mask_t, hRe)
        hIm = np.ma.masked_where(~srcgrd.mask_t, hIm)
        uRe = np.ma.masked_where(~srcgrd.mask_u, uRe)
        uIm = np.ma.masked_where(~srcgrd.mask_u, uIm)
        vRe = np.ma.masked_where(~srcgrd.mask_v, vRe)
        vIm = np.ma.masked_where(~srcgrd.mask_v, vIm)

        # remap h, u, v
        hRe = pyroms.remapping.remap(hRe, wts_file_t, spval=missing_value)
        hIm = pyroms.remapping.remap(hIm, wts_file_t, spval=missing_value)
        uRe = pyroms.remapping.remap(uRe, wts_file_u, spval=missing_value)
        uIm = pyroms.remapping.remap(uIm, wts_file_u, spval=missing_value)
        vRe = pyroms.remapping.remap(vRe, wts_file_v, spval=missing_value)
        vIm = pyroms.remapping.remap(vIm, wts_file_v, spval=missing_value)

    elif cst in consts2:
        # mask invalid data
        hRe = np.ma.masked_where(~srcgrd_lr.mask_t, hRe)
        hIm = np.ma.masked_where(~srcgrd_lr.mask_t, hIm)
        uRe = np.ma.masked_where(~srcgrd_lr.mask_u, uRe)
        uIm = np.ma.masked_where(~srcgrd_lr.mask_u, uIm)
        vRe = np.ma.masked_where(~srcgrd_lr.mask_v, vRe)
        vIm = np.ma.masked_where(~srcgrd_lr.mask_v, vIm)

        # remap h, u, v
        hRe = pyroms.remapping.remap(hRe, wts_file_t_lr, spval=missing_value)
        hIm = pyroms.remapping.remap(hIm, wts_file_t_lr, spval=missing_value)
        uRe = pyroms.remapping.remap(uRe, wts_file_u_lr, spval=missing_value)
        uIm = pyroms.remapping.remap(uIm, wts_file_u_lr, spval=missing_value)
        vRe = pyroms.remapping.remap(vRe, wts_file_v_lr, spval=missing_value)
        vIm = pyroms.remapping.remap(vIm, wts_file_v_lr, spval=missing_value)

    hamp[k, :, :] = (hRe**2+hIm**2)**0.5  # mm
    hpha[k, :, :] = np.arctan(-hIm/hRe)/np.pi*180.  # deg
    uamp = (uRe**2+uIm**2)**0.5  # mm
    upha = np.arctan(-uIm/uRe)/np.pi*180.  # deg
    vamp = (vRe**2+vIm**2)**0.5  # mm
    vpha = np.arctan(-vIm/vRe)/np.pi*180.  # deg


    # convert ap to ep
    cmax[k, :, :], ecc, cang[k, :, :], cpha[k, :, :], w = CGrid_TPXO8.tidal_ellipse.ap2ep(uamp, upha, vamp, vpha)
    cmin[k, :, :] = cmax[k, :, :]*ecc

# -------------------------------------------------------------------------
# mask land grid points
cmax[:, dstgrd.hgrid.mask_rho==0] = missing_value
cmin[:, dstgrd.hgrid.mask_rho==0] = missing_value
cang[:, dstgrd.hgrid.mask_rho==0] = missing_value
cpha[:, dstgrd.hgrid.mask_rho==0] = missing_value

# -------------------------------------------------------------------------
savedata = 1
if savedata == 1:
    # write tidal information to nc file
    # -------------------------------------------------------------------------
    # define tidal constituents names and periods
    tide_name = np.array([ list('Q1'), list('O1'), list('P1'), list('K1'),
                           list('N2'), list('M2'), list('S2'), list('K2'),
                           list('MF')])
    tide_period = np.array([26.8683567047119, 25.8193397521973, 24.0658893585205, 23.9344692230225,
                            12.6583499908447, 12.420599937439, 12, 11.9672346115112,
                            13.66079*24])

    # -------------------------------------------------------------------------
    # create nc file
    fh = nc.Dataset('Chukchi2_tides_otps.nc', 'w')
    fh.createDimension('namelen', 4)
    fh.createDimension('tide_period', consts_num)
    fh.createDimension('eta_rho', eta)
    fh.createDimension('xi_rho', xi)


    fh.history = 'Tides from TPXO8'
    import time
    fh.creation_date = time.strftime('%c')
    fh.Type = 'ROMS Tidal Forcing File'
    fh.Title = 'Forcing for' + dstgrd.name + 'domain'
    fh.grid_file = dstgrd.name
    fh.Source = 'OTPS'


    name_nc = fh.createVariable('tide_name', 'c', ('tide_period', 'namelen'))

    period_nc = fh.createVariable('tide_period', 'd', ('tide_period'), fill_value=-9999)
    period_nc.field = 'tide_period, scalar'
    period_nc.long_name = 'tidal angular period'
    period_nc.units = 'hours'

    lat_nc = fh.createVariable('lat_rho', 'd', ('eta_rho', 'xi_rho'))
    lat_nc.field = 'lat_rho, scalar'
    lat_nc.long_name = 'latitude of RHO-points'
    lat_nc.units = 'degree north'

    lon_nc = fh.createVariable('lon_rho', 'd', ('eta_rho', 'xi_rho'))
    lon_nc.field = 'lon_rho, scalar'
    lon_nc.long_name = 'longitude of RHO-points'
    lon_nc.units = 'degree east'

    msk_nc = fh.createVariable('mask_rho', 'd', ('eta_rho', 'xi_rho'))
    msk_nc.long_name = 'mask on RHO-points'
    msk_nc.option_0 = 'land'
    msk_nc.option_1 = 'water'

    Eamp_nc = fh.createVariable('tide_Eamp', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
    Eamp_nc.field = 'tide_Eamp, scalar'
    Eamp_nc.long_name = 'tidal elevation amplitude'
    Eamp_nc.units = 'meter'

    Ephase_nc = fh.createVariable('tide_Ephase', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
    Ephase_nc.field = 'tide_Ephase, scalar'
    Ephase_nc.long_name = 'tidal elevation phase angle'
    Ephase_nc.units = 'degrees, time of maximum elevation with respect chosen time orgin'

    Cmax_nc = fh.createVariable('tide_Cmax', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
    Cmax_nc.field = 'tide_Cmax, scalar'
    Cmax_nc.long_name = 'maximum tidal current, ellipse semi-major axis'
    Cmax_nc.units = 'meter second-1'

    Cmin_nc = fh.createVariable('tide_Cmin', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
    Cmin_nc.field = 'tide_Cmin, scalar'
    Cmin_nc.long_name = 'minimum tidal current, ellipse semi-minor axis'
    Cmin_nc.units = 'meter second-1'

    Cangle_nc = fh.createVariable('tide_Cangle', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
    Cangle_nc.field = 'tide_Cangle, scalar'
    Cangle_nc.long_name = 'tidal current inclination angle'
    Cangle_nc.units = 'degrees between semi-major axis and East'

    Cphase_nc = fh.createVariable('tide_Cphase', 'd', ('tide_period', 'eta_rho', 'xi_rho'), fill_value=-9999)
    Cphase_nc.field = 'tide_Cphase, scalar'
    Cphase_nc.long_name = 'tidal current phase angle'
    Cphase_nc.units = 'degrees, time of maximum velocity'

    # -------------------------------------------------------------------------
    # write data
    name_nc[:, 0:2] = tide_name
    period_nc[:] = tide_period

    lat_nc[:, :] = lat
    lon_nc[:, :] = lon
    msk_nc[:, :] = msk
    Eamp_nc[:, :, :] = hamp
    Ephase_nc[:, :, :] = hpha
    Cmax_nc[:, :, :] = cmax
    Cmin_nc[:, :, :] = cmin
    Cangle_nc[:, :, :] = cang
    Cphase_nc[:, :, :] = cpha

    fh.close()

