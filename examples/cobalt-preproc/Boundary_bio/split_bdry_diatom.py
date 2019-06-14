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
spval = -10000000000.

nlg = nc.variables['nlg_north'][:]
nlg *= 0.5
nc.variables['nlg_north'][:] = nlg

nc.createVariable('nmd_north', 'f8', ('ocean_time', 's_rho', 'xi_rho'), fill_value=spval)
nc.variables['nmd_north'].long_name = 'Medium Phytoplankton Nitrogen north boundary condition'
nc.variables['nmd_north'].units = 'mol/kg'
nc.variables['nmd_north'].time = 'ocean_time'
nc.variables['nmd_north'].field = 'nmd_north, scalar, series'
nc.variables['nmd_north'][:] = nlg

felg = nc.variables['felg_north'][:]
felg *= 0.5
nc.variables['felg_north'][:] = felg

nc.createVariable('femd_north', 'f8', ('ocean_time', 's_rho', 'xi_rho'), fill_value=spval)
nc.variables['femd_north'].long_name = 'Medium Phytoplankton Iron north boundary condition'
nc.variables['femd_north'].units = 'mol/kg'
nc.variables['femd_north'].time = 'ocean_time'
nc.variables['femd_north'].field = 'femd_north, scalar, series'
nc.variables['femd_north'][:] = felg

silg = nc.variables['silg_north'][:]
silg *= 0.5
nc.variables['silg_north'][:] = silg

nc.createVariable('simd_north', 'f8', ('ocean_time', 's_rho', 'xi_rho'), fill_value=spval)
nc.variables['simd_north'].long_name = 'Medium Phytoplankton Silicon north boundary condition'
nc.variables['simd_north'].units = 'mol/kg'
nc.variables['simd_north'].time = 'ocean_time'
nc.variables['simd_north'].field = 'simd_north, scalar, series'
nc.variables['simd_north'][:] = silg


nlg = nc.variables['nlg_south'][:]
nlg *= 0.5
nc.variables['nlg_south'][:] = nlg

nc.createVariable('nmd_south', 'f8', ('ocean_time', 's_rho', 'xi_rho'), fill_value=spval)
nc.variables['nmd_south'].long_name = 'Medium Phytoplankton Nitrogen south boundary condition'
nc.variables['nmd_south'].units = 'mol/kg'
nc.variables['nmd_south'].time = 'ocean_time'
nc.variables['nmd_south'].field = 'nmd_south, scalar, series'
nc.variables['nmd_south'][:] = nlg

felg = nc.variables['felg_south'][:]
felg *= 0.5
nc.variables['felg_south'][:] = felg

nc.createVariable('femd_south', 'f8', ('ocean_time', 's_rho', 'xi_rho'), fill_value=spval)
nc.variables['femd_south'].long_name = 'Medium Phytoplankton Iron south boundary condition'
nc.variables['femd_south'].units = 'mol/kg'
nc.variables['femd_south'].time = 'ocean_time'
nc.variables['femd_south'].field = 'femd_south, scalar, series'
nc.variables['femd_south'][:] = felg

silg = nc.variables['silg_south'][:]
silg *= 0.5
nc.variables['silg_south'][:] = silg

nc.createVariable('simd_south', 'f8', ('ocean_time', 's_rho', 'xi_rho'), fill_value=spval)
nc.variables['simd_south'].long_name = 'Medium Phytoplankton Silicon south boundary condition'
nc.variables['simd_south'].units = 'mol/kg'
nc.variables['simd_south'].time = 'ocean_time'
nc.variables['simd_south'].field = 'simd_south, scalar, series'
nc.variables['simd_south'][:] = silg

nlg = nc.variables['nlg_west'][:]
nlg *= 0.5
nc.variables['nlg_west'][:] = nlg

nc.createVariable('nmd_west', 'f8', ('ocean_time', 's_rho', 'eta_rho'), fill_value=spval)
nc.variables['nmd_west'].long_name = 'Medium Phytoplankton Nitrogen west boundary condition'
nc.variables['nmd_west'].units = 'mol/kg'
nc.variables['nmd_west'].time = 'ocean_time'
nc.variables['nmd_west'].field = 'nmd_west, scalar, series'
nc.variables['nmd_west'][:] = nlg

felg = nc.variables['felg_west'][:]
felg *= 0.5
nc.variables['felg_west'][:] = felg

nc.createVariable('femd_west', 'f8', ('ocean_time', 's_rho', 'eta_rho'), fill_value=spval)
nc.variables['femd_west'].long_name = 'Medium Phytoplankton Iron west boundary condition'
nc.variables['femd_west'].units = 'mol/kg'
nc.variables['femd_west'].time = 'ocean_time'
nc.variables['femd_west'].field = 'femd_west, scalar, series'
nc.variables['femd_west'][:] = felg

silg = nc.variables['silg_west'][:]
silg *= 0.5
nc.variables['silg_west'][:] = silg

nc.createVariable('simd_west', 'f8', ('ocean_time', 's_rho', 'eta_rho'), fill_value=spval)
nc.variables['simd_west'].long_name = 'Medium Phytoplankton Silicon west boundary condition'
nc.variables['simd_west'].units = 'mol/kg'
nc.variables['simd_west'].time = 'ocean_time'
nc.variables['simd_west'].field = 'simd_west, scalar, series'
nc.variables['simd_west'][:] = silg

nlg = nc.variables['nlg_east'][:]
nlg *= 0.5
nc.variables['nlg_east'][:] = nlg

nc.createVariable('nmd_east', 'f8', ('ocean_time', 's_rho', 'eta_rho'), fill_value=spval)
nc.variables['nmd_east'].long_name = 'Medium Phytoplankton Nitrogen east boundary condition'
nc.variables['nmd_east'].units = 'mol/kg'
nc.variables['nmd_east'].time = 'ocean_time'
nc.variables['nmd_east'].field = 'nmd_east, scalar, series'
nc.variables['nmd_east'][:] = nlg

felg = nc.variables['felg_east'][:]
felg *= 0.5
nc.variables['felg_east'][:] = felg

nc.createVariable('femd_east', 'f8', ('ocean_time', 's_rho', 'eta_rho'), fill_value=spval)
nc.variables['femd_east'].long_name = 'Medium Phytoplankton Iron east boundary condition'
nc.variables['femd_east'].units = 'mol/kg'
nc.variables['femd_east'].time = 'ocean_time'
nc.variables['femd_east'].field = 'femd_east, scalar, series'
nc.variables['femd_east'][:] = felg

silg = nc.variables['silg_east'][:]
silg *= 0.5
nc.variables['silg_east'][:] = silg

nc.createVariable('simd_east', 'f8', ('ocean_time', 's_rho', 'eta_rho'), fill_value=spval)
nc.variables['simd_east'].long_name = 'Medium Phytoplankton Silicon east boundary condition'
nc.variables['simd_east'].units = 'mol/kg'
nc.variables['simd_east'].time = 'ocean_time'
nc.variables['simd_east'].field = 'simd_east, scalar, series'
nc.variables['simd_east'][:] = silg

nc.close()
