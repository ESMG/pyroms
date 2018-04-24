import subprocess
import os
import sys
import subprocess
import numpy as np

import pyroms
import pyroms_toolbox

from remap_bio_woa import remap_bio_woa
from remap_bio_glodap import remap_bio_glodap

data_dir_woa = '/archive/AKWATERS/kshedstrom/COBALT/'
data_dir_glodap = '/archive/AKWATERS/kshedstrom/COBALT/'
dst_dir='./'

src_grd = pyroms_toolbox.BGrid_GFDL.get_nc_BGrid_GFDL('/archive/AKWATERS/kshedstrom/COBALT/GFDL_CM2.1_grid.nc', \
          name='ESM2M_ARCTIC4', area='tripole', ystart=150)
dst_grd = pyroms.grid.get_ROMS_grid('ARCTIC4')

# tracer informations
list_tracer = ['alk', 'cadet_arag', 'cadet_calc', 'dic', 'fed', 'fedet', 'fedi', 'felg', 'fesm', 'ldon', 'ldop', 'lith', 'lithdet', 'nbact', 'ndet', 'ndi', 'nlg', 'nsm', 'nh4', 'no3', 'o2', 'pdet', 'po4', 'srdon', 'srdop', 'sldon', 'sldop', 'sidet', 'silg', 'sio4', 'nsmz', 'nmdz', 'nlgz','cased','chl','irr_mem','htotal','co3_ion']

tracer_longname = ['Alkalinity', 'Detrital CaCO3', 'Detrital CaCO3', 'Dissolved Inorganic Carbon', 'Dissolved Iron', 'Detrital Iron', 'Diazotroph Iron', 'Large Phytoplankton Iron', 'Small Phytoplankton Iron', 'labile DON', 'labile DOP', 'Lithogenic Aluminosilicate', 'lithdet', 'bacterial', 'ndet', 'Diazotroph Nitrogen', 'Large Phytoplankton Nitrogen', 'Small Phytoplankton Nitrogen', 'Ammonia', 'Nitrate', 'Oxygen', 'Detrital Phosphorus', 'Phosphate', 'Semi-Refractory DON', 'Semi-Refractory DOP', 'Semilabile DON', 'Semilabile DOP', 'Detrital Silicon', 'Large Phytoplankton Silicon', 'Silicate', 'Small Zooplankton Nitrogen', 'Medium-sized zooplankton Nitrogen', 'large Zooplankton Nitrogen','Sediment CaCO3','Cholorophyll','Irradiance Memory','Total H+','Carbonate ion']

tracer_units = ['mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'g/kg', 'g/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg', 'mol/kg','mol.m-3','ug.kg-1','W.m-2','mol/kg','mol/kg']


ic_file = dst_dir + dst_grd.name + '_ic_bio_GFDL-APR.nc'

#------- WOA13 ---------------------------------
id_tracer_update_woa = [19,20,22,29]
list_tracer_update_woa = []
tracer_longname_update_woa = []
tracer_units_update_woa = []

for idtra in id_tracer_update_woa:
        print(list_tracer[idtra])

for idtra in id_tracer_update_woa:
        # add to tracer update
        list_tracer_update_woa.append(list_tracer[idtra])
        tracer_longname_update_woa.append(tracer_longname[idtra])
        tracer_units_update_woa.append(tracer_units[idtra])

for ktr in np.arange(len(list_tracer_update_woa)):
    ctra = list_tracer_update_woa[ktr]
    if ctra == 'sio4':
       ctra = 'si'
    mydict = {'tracer':list_tracer_update_woa[ktr],'longname':tracer_longname_update_woa[ktr],'units':tracer_units_update_woa[ktr],'file':data_dir_woa + ctra + '_WOA13-CM2.1_monthly.nc', \
    'nframe':3}
    remap_bio_woa(mydict, src_grd, dst_grd, dst_dir=dst_dir)
    out_file = dst_dir + dst_grd.name + '_ic_bio_' + list_tracer_update_woa[ktr] + '.nc'
    command = ('ncks', '-a', '-A', out_file, ic_file)
    subprocess.check_call(command)
    os.remove(out_file)

#--------- GLODAP -------------------------------
id_tracer_update_glodap = [0,3]
list_tracer_update_glodap = []
tracer_longname_update_glodap = []
tracer_units_update_glodap = []

for idtra in id_tracer_update_glodap:
        print(list_tracer[idtra])

for idtra in id_tracer_update_glodap:
        # add to tracer update
        list_tracer_update_glodap.append(list_tracer[idtra])
        tracer_longname_update_glodap.append(tracer_longname[idtra])
        tracer_units_update_glodap.append(tracer_units[idtra])

for ktr in np.arange(len(list_tracer_update_glodap)):
    ctra = list_tracer_update_glodap[ktr]
    mydict = {'tracer':list_tracer_update_glodap[ktr],'longname':tracer_longname_update_glodap[ktr],'units':tracer_units_update_glodap[ktr],'file':data_dir_glodap + ctra + '_GLODAP-ESM2M_annual.nc', \
    'nframe':0}
    remap_bio_glodap(mydict, src_grd, dst_grd, dst_dir=dst_dir)
    out_file = dst_dir + dst_grd.name + '_ic_bio_' + list_tracer_update_glodap[ktr] + '.nc'
    command = ('ncks', '-a', '-A', out_file, ic_file)
    subprocess.check_call(command)
    os.remove(out_file)


