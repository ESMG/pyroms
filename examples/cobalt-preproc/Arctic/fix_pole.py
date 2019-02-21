import subprocess
import os
import sys
import subprocess
import numpy as np
import netCDF4 as nc

dst_dir='./'

ic_file = dst_dir + 'ARCTIC4_ic_bio_GFDL-APR.nc'
fidic = nc.Dataset(ic_file,'a')
Cs_r = fidic.variables['Cs_r']
nz = Cs_r.shape[0]

# define all tracer stuff
list_tracer = ['alk', 'cadet_arag', 'cadet_calc', 'dic', 'fed', 'fedet', 'fedi', 'felg', 'fesm', 'ldon', 'ldop', 'lith', 'lithdet', 'nbact', 'ndet', 'ndi', 'nlg', 'nsm', 'nh4', 'no3', 'o2', 'pdet', 'po4', 'srdon', 'srdop', 'sldon', 'sldop', 'sidet', 'silg', 'sio4', 'nsmz', 'nmdz', 'nlgz','cased','chl','irr_mem','htotal','co3_ion']


print('\nFixing a north pole problem')
for tr in list_tracer:
    print('for variable', tr)
    tracer = fidic.variables[tr][:]
    mysum = np.zeros((nz))
    count = 0
    for j in range(753,768):
        for i in range(271,287):
            if tracer[0,0,j,i] != 0:
                count += 1
                mysum += tracer[0,:,j,i]
    print('count', count)
    mysum = mysum/count
    print('mysum', mysum)
    for j in range(753,768):
        for i in range(271,287):
            if tracer[0,0,j,i] == 0:
                tracer[0,:,j,i] = mysum
    fidic.variables[tr][:] = tracer

# These two tracers contain zeros, leading to nans.
tracer = fidic.variables['cased'][:]
mysum = 0.25*(tracer[0,:,752,279] + tracer[0,:,768,279] + tracer[0,:,760,270] + tracer[0,:,602,287])
for j in range(753,768):
    for i in range(271,287):
        tracer[0,:,j,i] = mysum
fidic.variables['cased'][:] = tracer

tracer = fidic.variables['irr_mem'][:]
mysum = 0.25*(tracer[0,:,752,279] + tracer[0,:,768,279] + tracer[0,:,760,270] + tracer[0,:,602,287])
for j in range(753,768):
    for i in range(271,287):
        tracer[0,:,j,i] = mysum
fidic.variables['irr_mem'][:] = tracer

fidic.close()
