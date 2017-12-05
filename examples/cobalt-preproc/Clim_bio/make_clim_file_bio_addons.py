import subprocess
import os
import sys
import commands
import numpy as np

import pyroms
import pyroms_toolbox

from remap_bio_woa import remap_bio_woa
from remap_bio_glodap import remap_bio_glodap

data_dir_woa = '/archive/u1/uaf/kate/COBALT/'
data_dir_glodap = '/archive/u1/uaf/kate/COBALT/'
dst_dir='./'

src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL('/archive/u1/uaf/kate/COBALT/GFDL_CM2.1_grid.nc', name='ESM2M_NWGOA3')
dst_grd = pyroms.grid.get_ROMS_grid('NWGOA3')

# define all tracer stuff
list_tracer = ['alk', 'cadet_arag', 'cadet_calc', 'dic', 'fed', 'fedet', 'fedi', 'felg', 'fesm', 'ldon', 'ldop', 'lith', 'lithdet', 'nbact', 'ndet', 'ndi', 'nlg', 'nsm', 'nh4', 'no3', 'o2', 'pdet', 'po4', 'srdon', 'srdop', 'sldon', 'sldop', 'sidet', 'silg', 'sio4', 'nsmz', 'nmdz', 'nlgz']

tracer_longname = ['Alkalinity', 'Detrital CaCO3', 'Detrital CaCO3', 'Dissolved Inorganic Carbon', 'Dissolved Iron', 'Detrital Iron', 'Diazotroph Iron', 'Large Phytoplankton Iron', 'Small Phytoplankton Iron', 'labile DON', 'labile DOP', 'Lithogenic Aluminosilicate', 'lithdet', 'bacterial', 'ndet', 'Diazotroph Nitrogen', 'Large Phytoplankton Nitrogen', 'Small Phytoplankton Nitrogen', 'Ammonia', 'Nitrate', 'Oxygen', 'Detrital Phosphorus', 'Phosphate', 'Semi-Refractory DON', 'Semi-Refractory DOP', 'Semilabile DON', 'Semilabile DOP', 'Detrital Silicon', 'Large Phytoplankton Silicon', 'Silicate', 'Small Zooplankton Nitrogen', 'Medium-sized zooplankton Nitrogen', 'large Zooplankton Nitrogen']

tracer_units = ['mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'g/kg', 'g/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg']


#------- WOA13 ---------------------------------
id_tracer_update_woa = [19,20,22,29]
list_tracer_update_woa = []
tracer_longname_update_woa = []
tracer_units_update_woa = []

for idtra in id_tracer_update_woa:
        print list_tracer[idtra]

for idtra in id_tracer_update_woa:
        # add to tracer update
        list_tracer_update_woa.append(list_tracer[idtra])
        tracer_longname_update_woa.append(tracer_longname[idtra])
        tracer_units_update_woa.append(tracer_units[idtra])

for mm in np.arange(12):
    clim_file = dst_dir + dst_grd.name + '_clim_bio_GFDL+WOA+GLODAP_m' + str(mm+1).zfill(2) + '.nc'
    print '\nBuild CLIM file for month', mm
    for ktr in np.arange(len(list_tracer_update_woa)):
        ctra = list_tracer_update_woa[ktr]
        if ctra == 'sio4':
           ctra = 'si'
        mydict = {'tracer':list_tracer_update_woa[ktr],'longname':tracer_longname_update_woa[ktr],'units':tracer_units_update_woa[ktr],'file':data_dir_woa + ctra + '_WOA13-CM2.1_monthly.nc', \
        'frame':mm}
        remap_bio_woa(mydict, src_grd, dst_grd, dst_dir=dst_dir)
        out_file = dst_dir + dst_grd.name + '_clim_bio_' + list_tracer_update_woa[ktr] + '.nc'
        command = ('ncks', '-a', '-A', out_file, clim_file)
        subprocess.check_call(command)
        os.remove(out_file)

#--------- GLODAP -------------------------------
id_tracer_update_glodap = [0,3]
list_tracer_update_glodap = []
tracer_longname_update_glodap = []
tracer_units_update_glodap = []

for idtra in id_tracer_update_glodap:
        print list_tracer[idtra]

for idtra in id_tracer_update_glodap:
        # add to tracer update
        list_tracer_update_glodap.append(list_tracer[idtra])
        tracer_longname_update_glodap.append(tracer_longname[idtra])
        tracer_units_update_glodap.append(tracer_units[idtra])

for mm in np.arange(12):
    clim_file = dst_dir + dst_grd.name + '_clim_bio_GFDL+WOA+GLODAP_m' + str(mm+1).zfill(2) + '.nc'
    print '\nBuild CLIM file for month', mm
    for ktr in np.arange(len(list_tracer_update_glodap)):
        ctra = list_tracer_update_glodap[ktr]
        mydict = {'tracer':list_tracer_update_glodap[ktr],'longname':tracer_longname_update_glodap[ktr],'units':tracer_units_update_glodap[ktr],'file':data_dir_glodap + ctra + '_GLODAP-ESM2M_annual.nc', \
        'frame':mm}
        remap_bio_glodap(mydict, src_grd, dst_grd, dst_dir=dst_dir)
        out_file = dst_dir + dst_grd.name + '_clim_bio_' + list_tracer_update_glodap[ktr] + '.nc'
        command = ('ncks', '-a', '-A', out_file, clim_file)
        subprocess.check_call(command)
        os.remove(out_file)

