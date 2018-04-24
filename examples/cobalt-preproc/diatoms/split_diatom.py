import numpy as np
import netCDF4
import sys

# This program splits the large phytoplankton into large and medium.
# True for nitrogen, iron, and silicon.
# Copy the original file before running this, then operate on the copy.
# Could probably do this with ncap2 instead.
#
ncfile = sys.argv[1]
nc = netCDF4.Dataset(ncfile, 'a', format='NETCDF3_CLASSIC')
spval = 1.e+37

nlg = nc.variables['nlg'][:]
nlg *= 0.5
nc.variables['nlg'][:] = nlg

nc.createVariable('nmd', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
#nc.createVariable('nmd', 'f8', ('ocean_time', 'two', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
nc.variables['nmd'].long_name = 'Medium Phytoplankton Nitrogen'
nc.variables['nmd'].units = 'mol/kg'
nc.variables['nmd'].time = 'ocean_time'
nc.variables['nmd'].field = 'nmd, scalar, series'
nc.variables['nmd'].grid = 'grid'
nc.variables['nmd'].location = 'face'
nc.variables['nmd'].coordinates = 'lon_rho lat_rho s_rho ocean_time'
nc.variables['nmd'][:] = nlg

felg = nc.variables['felg'][:]
felg *= 0.5
nc.variables['felg'][:] = felg

nc.createVariable('femd', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
#nc.createVariable('femd', 'f8', ('ocean_time', 'two', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
nc.variables['femd'].long_name = 'Medium Phytoplankton Iron'
nc.variables['femd'].units = 'mol/kg'
nc.variables['femd'].time = 'ocean_time'
nc.variables['femd'].field = 'femd, scalar, series'
nc.variables['femd'].grid = 'grid'
nc.variables['femd'].location = 'face'
nc.variables['femd'].coordinates = 'lon_rho lat_rho s_rho ocean_time'
nc.variables['femd'][:] = felg

silg = nc.variables['silg'][:]
silg *= 0.5
nc.variables['silg'][:] = silg

nc.createVariable('simd', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
#nc.createVariable('simd', 'f8', ('ocean_time', 'two', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
nc.variables['simd'].long_name = 'Medium Phytoplankton Silicon'
nc.variables['simd'].units = 'mol/kg'
nc.variables['simd'].time = 'ocean_time'
nc.variables['simd'].field = 'simd, scalar, series'
nc.variables['simd'].grid = 'grid'
nc.variables['simd'].location = 'face'
nc.variables['simd'].coordinates = 'lon_rho lat_rho s_rho ocean_time'
nc.variables['simd'][:] = silg

mu_mem_lg = nc.variables['mu_mem_lg'][:]
mu_mem_lg *= 0.5
nc.variables['mu_mem_lg'][:] = mu_mem_lg

nc.createVariable('mu_mem_md', 'f8', ('two', 's_rho', 'eta_rho', 'xi_rho'), fill_value=spval)
nc.variables['mu_mem_md'].long_name = 'Medium Phytoplankton Silicon'
nc.variables['mu_mem_md'].units = 'mol/kg'
nc.variables['mu_mem_md'].time = 'ocean_time'
nc.variables['mu_mem_md'].field = 'mu_mem_md, scalar, series'
nc.variables['mu_mem_md'].grid = 'grid'
nc.variables['mu_mem_md'].location = 'face'
nc.variables['mu_mem_md'].coordinates = 'lon_rho lat_rho s_rho ocean_time'
nc.variables['mu_mem_md'][:] = mu_mem_lg

nc.close()
