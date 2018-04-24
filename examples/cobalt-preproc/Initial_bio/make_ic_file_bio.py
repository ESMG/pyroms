import subprocess
import os
import sys
import subprocess
import numpy as np
import netCDF4 as nc

import pyroms
import pyroms_toolbox

from remap_bio import remap_bio

#build list of date to remap
tag = 'y1988-2007m01'

data_dir = '/archive/u1/uaf/kate/COBALT/'
dst_dir='./'

src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL('/archive/u1/uaf/kate/COBALT/GFDL_CM2.1_grid.nc', name='ESM2M_NWGOA3')
dst_grd = pyroms.grid.get_ROMS_grid('NWGOA3')

# define all tracer stuff
list_tracer = ['alk', 'cadet_arag', 'cadet_calc', 'dic', 'fed', 'fedet', 'fedi', 'felg', 'fesm', 'ldon', 'ldop', 'lith', 'lithdet', 'nbact', 'ndet', 'ndi', 'nlg', 'nsm', 'nh4', 'no3', 'o2', 'pdet', 'po4', 'srdon', 'srdop', 'sldon', 'sldop', 'sidet', 'silg', 'sio4', 'nsmz', 'nmdz', 'nlgz','cased','chl','irr_mem','htotal','co3_ion']

tracer_longname = ['Alkalinity', 'Detrital CaCO3', 'Detrital CaCO3', 'Dissolved Inorganic Carbon', 'Dissolved Iron', 'Detrital Iron', 'Diazotroph Iron', 'Large Phytoplankton Iron', 'Small Phytoplankton Iron', 'labile DON', 'labile DOP', 'Lithogenic Aluminosilicate', 'lithdet', 'bacterial', 'ndet', 'Diazotroph Nitrogen', 'Large Phytoplankton Nitrogen', 'Small Phytoplankton Nitrogen', 'Ammonia', 'Nitrate', 'Oxygen', 'Detrital Phosphorus', 'Phosphate', 'Semi-Refractory DON', 'Semi-Refractory DOP', 'Semilabile DON', 'Semilabile DOP', 'Detrital Silicon', 'Large Phytoplankton Silicon', 'Silicate', 'Small Zooplankton Nitrogen', 'Medium-sized zooplankton Nitrogen', 'large Zooplankton Nitrogen','Sediment CaCO3','Cholorophyll','Irradiance Memory','Total H+','Carbonate ion']

tracer_units = ['mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'g/kg', 'g/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg','mol.m-3','ug.kg-1','W.m-2','mol/kg','mol/kg']



print('\nBuild IC file for time %s' %tag)
for ktr in np.arange(len(list_tracer)):
    mydict = {'tracer':list_tracer[ktr],'longname':tracer_longname[ktr],'units':tracer_units[ktr], \
        'file':data_dir + 'ocean_cobalt_tracers.1988-2007.01_12.nc', 'nframe':0}
    remap_bio(mydict, src_grd, dst_grd, dst_dir=dst_dir)

## merge file
ic_file = dst_dir + dst_grd.name + '_ic_bio_GFDL-JAN.nc'
out_file = dst_dir + dst_grd.name + '_ic_bio_' + list_tracer[0] + '.nc'
command = ('ncks', '-a', '-O', out_file, ic_file)
subprocess.check_call(command)
os.remove(out_file)

for ktr in np.arange(1,len(list_tracer)):
    out_file = dst_dir + dst_grd.name + '_ic_bio_' + list_tracer[ktr] + '.nc'
    command = ('ncks', '-a', '-A', out_file, ic_file)
    subprocess.check_call(command)
    os.remove(out_file)


#------------------ Add additional zeros fields ----------------------------------

fidic = nc.Dataset(ic_file,'a',format='NETCDF3_64BIT')

fidic.createVariable('mu_mem_lg', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',) )
fidic.variables['mu_mem_lg'].long_name = 'large phytoplankton aggregation memory'
fidic.variables['mu_mem_lg'].units = ''
fidic.variables['mu_mem_lg'].field = 'mu_mem_lg, scalar, series'

fidic.createVariable('mu_mem_di', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',) )
fidic.variables['mu_mem_di'].long_name = 'medium phytoplankton aggregation memory'
fidic.variables['mu_mem_di'].units = ''
fidic.variables['mu_mem_di'].field = 'mu_mem_di, scalar, series'

fidic.createVariable('mu_mem_sm', 'f8', ('ocean_time', 's_rho', 'eta_rho', 'xi_rho',) )
fidic.variables['mu_mem_sm'].long_name = 'small phytoplankton aggregation memory'
fidic.variables['mu_mem_sm'].units = ''
fidic.variables['mu_mem_sm'].field = 'mu_mem_sm, scalar, series'

fidic.variables['mu_mem_lg'][0,:,:,:] = 0.
fidic.variables['mu_mem_di'][0,:,:,:] = 0.
fidic.variables['mu_mem_sm'][0,:,:,:] = 0.

fidic.close()
