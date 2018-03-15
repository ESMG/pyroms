import numpy as np
from datetime import datetime
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF


def nc_create_roms_file(filename, grd, ocean_time, lgrid=True):

    # create file
    nc = netCDF.Dataset(filename, 'w', format='NETCDF3_64BIT')
    nc.Description = 'ROMS file'
    nc.Author = 'pyroms_toolbox.nc_create_roms_file'
    nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    nc.title = 'ROMS file'

    nc.createDimension('xi_rho', np.size(grd.hgrid.mask_rho,1))
    nc.createDimension('xi_u', np.size(grd.hgrid.mask_u,1))
    nc.createDimension('xi_v', np.size(grd.hgrid.mask_v,1))
    nc.createDimension('xi_psi', np.size(grd.hgrid.mask_psi,1))
    nc.createDimension('eta_rho', np.size(grd.hgrid.mask_rho,0))
    nc.createDimension('eta_u', np.size(grd.hgrid.mask_u,0))
    nc.createDimension('eta_v', np.size(grd.hgrid.mask_v,0))
    nc.createDimension('eta_psi', np.size(grd.hgrid.mask_psi,0))
    nc.createDimension('s_rho', grd.vgrid.N)
    nc.createDimension('s_w', grd.vgrid.Np)
    nc.createDimension('ocean_time', None)

    # write time and grid information
    nc.createVariable('theta_s', 'f8', ())
    nc.variables['theta_s'].long_name = 'S-coordinate surface control parameter'
    nc.variables['theta_s'][:] = grd.vgrid.theta_s

    nc.createVariable('theta_b', 'f8', ())
    nc.variables['theta_b'].long_name = 'S-coordinate bottom control parameter'
    nc.variables['theta_b'][:] = grd.vgrid.theta_b

    nc.createVariable('Tcline', 'f8', ())
    nc.variables['Tcline'].long_name = 'S-cordinate surface/bottom layer width'
    nc.variables['Tcline'].units = 'meter'
    nc.variables['Tcline'][:] = grd.vgrid.Tcline

    nc.createVariable('hc', 'f8', ())
    nc.variables['hc'].long_name = 'S-coordinate parameter, critical depth'
    nc.variables['hc'].units = 'meter'
    nc.variables['hc'][:] = grd.vgrid.hc

    nc.createVariable('s_rho', 'f8', ('s_rho'))
    nc.variables['s_rho'].long_name = 'S-coordinate at RHO-points'
    nc.variables['s_rho'].valid_min = '-1.0'
    nc.variables['s_rho'].valid_max = '0.0'
    nc.variables['s_rho'].field = 's_rho,scalar'
    nc.variables['s_rho'][:] = grd.vgrid.s_rho

    nc.createVariable('s_w', 'f8', ('s_w'))
    nc.variables['s_w'].long_name = 'S-coordinate at W-points'
    nc.variables['s_w'].valid_min = '-1.0'
    nc.variables['s_w'].valid_max = '0.0'
    nc.variables['s_w'].field = 's_w,scalar'
    nc.variables['s_w'][:] = grd.vgrid.s_w

    nc.createVariable('Cs_r', 'f8', ('s_rho'))
    nc.variables['Cs_r'].long_name = 'S-coordinate stretching curves at RHO-points'
    nc.variables['Cs_r'].valid_min = '-1.0'
    nc.variables['Cs_r'].valid_max = '0.0'
    nc.variables['Cs_r'].field = 'Cs_r,scalar'
    nc.variables['Cs_r'][:] = grd.vgrid.Cs_r

    nc.createVariable('Cs_w', 'f8', ('s_w'))
    nc.variables['Cs_w'].long_name = 'S-coordinate stretching curves at W-points'
    nc.variables['Cs_w'].valid_min = '-1.0'
    nc.variables['Cs_w'].valid_max = '0.0'
    nc.variables['Cs_w'].field = 'Cs_w,scalar'
    nc.variables['Cs_w'][:] = grd.vgrid.Cs_w

    if (lgrid):
        nc.createVariable('h', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['h'].long_name = 'bathymetry at RHO-points'
        nc.variables['h'].units ='meter'
        nc.variables['h'].coordinates = 'lon_rho lat_rho'
        nc.variables['h'].field = 'bath, scalar'
        nc.variables['h'][:] = grd.vgrid.h

        nc.createVariable('pm', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['pm'].long_name = 'curvilinear coordinate metric in XI'
        nc.variables['pm'].units ='meter-1'
        nc.variables['pm'].coordinates = 'lon_rho lat_rho'
        nc.variables['pm'].field = 'pm, scalar'
        nc.variables['pm'][:] = 1. / grd.hgrid.dx

        nc.createVariable('pn', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['pn'].long_name = 'curvilinear coordinate metric in ETA'
        nc.variables['pn'].units ='meter-1'
        nc.variables['pn'].coordinates = 'lon_rho lat_rho'
        nc.variables['pn'].field = 'pn, scalar'
        nc.variables['pn'][:] = 1. / grd.hgrid.dy

        nc.createVariable('lon_rho', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['lon_rho'].long_name = 'longitude of RHO-points'
        nc.variables['lon_rho'].units = 'degree_east'
        nc.variables['lon_rho'].field = 'lon_rho, scalar'
        nc.variables['lon_rho'][:] = grd.hgrid.lon_rho

        nc.createVariable('lat_rho', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['lat_rho'].long_name = 'latitude of RHO-points'
        nc.variables['lat_rho'].units = 'degree_north'
        nc.variables['lat_rho'].field = 'lat_rho, scalar'
        nc.variables['lat_rho'][:] = grd.hgrid.lat_rho

        nc.createVariable('lon_u', 'f8', ('eta_u', 'xi_u'))
        nc.variables['lon_u'].long_name = 'longitude of U-points'
        nc.variables['lon_u'].units = 'degree_east'
        nc.variables['lon_u'].field = 'lon_u, scalar'
        nc.variables['lon_u'][:] = grd.hgrid.lon_u

        nc.createVariable('lat_u', 'f8', ('eta_u', 'xi_u'))
        nc.variables['lat_u'].long_name = 'latitude of U-points'
        nc.variables['lat_u'].units = 'degree_north'
        nc.variables['lat_u'].field = 'lat_u, scalar'
        nc.variables['lat_u'][:] = grd.hgrid.lat_u

        nc.createVariable('lon_v', 'f8', ('eta_v', 'xi_v'))
        nc.variables['lon_v'].long_name = 'longitude of V-points'
        nc.variables['lon_v'].units = 'degree_east'
        nc.variables['lon_v'].field = 'lon_v, scalar'
        nc.variables['lon_v'][:] = grd.hgrid.lon_v

        nc.createVariable('lat_v', 'f8', ('eta_v', 'xi_v'))
        nc.variables['lat_v'].long_name = 'latitude of V-points'
        nc.variables['lat_v'].units = 'degree_north'
        nc.variables['lat_v'].field = 'lat_v, scalar'
        nc.variables['lat_v'][:] = grd.hgrid.lat_v

        nc.createVariable('lon_psi', 'f8', ('eta_psi', 'xi_psi'))
        nc.variables['lon_psi'].long_name = 'longitude of PSI-points'
        nc.variables['lon_psi'].units = 'degree_east'
        nc.variables['lon_psi'].field = 'lon_psi, scalar'
        nc.variables['lon_psi'][:] = grd.hgrid.lon_psi

        nc.createVariable('lat_psi', 'f8', ('eta_psi', 'xi_psi'))
        nc.variables['lat_psi'].long_name = 'latitude of PSI-points'
        nc.variables['lat_psi'].units = 'degree_north'
        nc.variables['lat_psi'].field = 'lat_psi, scalar'
        nc.variables['lat_psi'][:] = grd.hgrid.lat_psi

        nc.createVariable('angle', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['angle'].long_name = 'angle between XI-axis and EAST'
        nc.variables['angle'].units = 'radians'
        nc.variables['angle'].coordinates = 'lon_rho lat_rho'
        nc.variables['angle'].field = 'angle, scalar'
        nc.variables['angle'][:] = grd.hgrid.angle_rho

        nc.createVariable('mask_rho', 'f8', ('eta_rho', 'xi_rho'))
        nc.variables['mask_rho'].long_name = 'mask on RHO-points'
        nc.variables['mask_rho'].option_0 = 'land'
        nc.variables['mask_rho'].option_1 = 'water'
        nc.variables['mask_rho'].coordinates = 'lon_rho lat_rho'
        nc.variables['mask_rho'][:] = grd.hgrid.mask_rho

        nc.createVariable('mask_u', 'f8', ('eta_u', 'xi_u'))
        nc.variables['mask_u'].long_name = 'mask on U-points'
        nc.variables['mask_u'].option_0 = 'land'
        nc.variables['mask_u'].option_1 = 'water'
        nc.variables['mask_u'].coordinates = 'lon_u lat_u'
        nc.variables['mask_u'][:] = grd.hgrid.mask_u

        nc.createVariable('mask_v', 'f8', ('eta_v', 'xi_v'))
        nc.variables['mask_v'].long_name = 'mask on V-points'
        nc.variables['mask_v'].option_0 = 'land'
        nc.variables['mask_v'].option_1 = 'water'
        nc.variables['mask_v'].coordinates = 'lon_v lat_v'
        nc.variables['mask_v'][:] = grd.hgrid.mask_v

        nc.createVariable('mask_psi', 'f8', ('eta_psi', 'xi_psi'))
        nc.variables['mask_psi'].long_name = 'mask on PSI-points'
        nc.variables['mask_psi'].option_0 = 'land'
        nc.variables['mask_psi'].option_1 = 'water'
        nc.variables['mask_psi'].coordinates = 'lon_psi lat_psi'
        nc.variables['mask_psi'][:] = grd.hgrid.mask_psi

    nc.createVariable('ocean_time', 'f8', ('ocean_time'))
    nc.variables['ocean_time'].long_name = ocean_time.long_name
    nc.variables['ocean_time'].units = ocean_time.units
    try:
        nc.variables['ocean_time'].field = ocean_time.field
    except:
        nc.variables['ocean_time'].field = ' '

    nc.close()
