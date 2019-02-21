import subprocess
import os
import sys
import subprocess
import numpy as np

import pyroms
import pyroms_toolbox

from remap_bio import remap_bio

data_dir = '/archive/u1/uaf/kate/COBALT/'
dst_dir='./'

src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL('/archive/u1/uaf/kate/COBALT/GFDL_CM2.1_grid.nc', name='ESM2M_NWGOA3')
dst_grd = pyroms.grid.get_ROMS_grid('NWGOA3')

# define all tracer stuff
list_tracer = ['alk', 'cadet_arag', 'cadet_calc', 'dic', 'fed', 'fedet', 'fedi', 'felg', 'fesm', 'ldon', 'ldop', 'lith', 'lithdet', 'nbact', 'ndet', 'ndi', 'nlg', 'nsm', 'nh4', 'no3', 'o2', 'pdet', 'po4', 'srdon', 'srdop', 'sldon', 'sldop', 'sidet', 'silg', 'sio4', 'nsmz', 'nmdz', 'nlgz']

tracer_longname = ['Alkalinity', 'Detrital CaCO3', 'Detrital CaCO3', 'Dissolved Inorganic Carbon', 'Dissolved Iron', 'Detrital Iron', 'Diazotroph Iron', 'Large Phytoplankton Iron', 'Small Phytoplankton Iron', 'labile DON', 'labile DOP', 'Lithogenic Aluminosilicate', 'lithdet', 'bacterial', 'ndet', 'Diazotroph Nitrogen', 'Large Phytoplankton Nitrogen', 'Small Phytoplankton Nitrogen', 'Ammonia', 'Nitrate', 'Oxygen', 'Detrital Phosphorus', 'Phosphate', 'Semi-Refractory DON', 'Semi-Refractory DOP', 'Semilabile DON', 'Semilabile DOP', 'Detrital Silicon', 'Large Phytoplankton Silicon', 'Silicate', 'Small Zooplankton Nitrogen', 'Medium-sized zooplankton Nitrogen', 'large Zooplankton Nitrogen']

tracer_units = ['mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'g/kg', 'g/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg']

for mm in np.arange(12):
    print('\nBuild clim file for month ', mm)
    for ktr in np.arange(len(list_tracer)):
        mydict = {'tracer':list_tracer[ktr],'longname':tracer_longname[ktr],'units':tracer_units[ktr],'file':data_dir + 'ocean_cobalt_tracers.1988-2007.01_12.nc','frame':mm}
        remap_bio(mydict, src_grd, dst_grd, dst_dir=dst_dir)

    # merge file
    clim_file = dst_dir + dst_grd.name + '_clim_bio_GFDL+WOA+GLODAP_m' + str(mm+1).zfill(2) + '.nc' 
    out_file = dst_dir + dst_grd.name + '_clim_bio_' + list_tracer[0] + '.nc' 
    command = ('ncks', '-a', '-O', out_file, clim_file)
    subprocess.check_call(command)
    os.remove(out_file)

    for ktr in np.arange(1,len(list_tracer)):
        out_file = dst_dir + dst_grd.name + '_clim_bio_' + list_tracer[ktr] + '.nc' 
        command = ('ncks', '-a', '-A', out_file, clim_file)
        subprocess.check_call(command)
        os.remove(out_file)
