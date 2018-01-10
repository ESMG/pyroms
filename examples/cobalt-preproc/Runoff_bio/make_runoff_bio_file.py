import subprocess
import os
import sys
import subprocess
import numpy as np
from datetime import *

import pyroms
import pyroms_toolbox
import netCDF4 as nc


file_nutrients='/archive/u1/uaf/kate/NWGOA3/Files_cobalt/NEWS-NWGOA3_clim.nc'
file_runoff='/archive/u1/uaf/kate/NWGOA3/Files_cobalt/runoff_NWGOA3_daitren_clim.nc'

dst_dir='./'

grd = pyroms.grid.get_ROMS_grid('NWGOA3')

fid_r = nc.Dataset(file_runoff,'r')
fid_n = nc.Dataset(file_nutrients,'r')

# get time
nctime = fid_r.variables['runoff_time'][:]
nt = len(nctime)

# define all tracer stuff
list_tracer = ['NO3','LDON','SLDON','SRDON','NDET','PO4','LDOP','SLDOP','SRDOP']
tracer_longname = []
tracer_units = []
for kt in list_tracer:
    tracer_longname.append(kt.lower() + ' river source')
    tracer_units.append('mol.m-2.s-1')

for ktr in np.arange(len(list_tracer)):
    dst_varname = 'river_' + list_tracer[ktr].lower()
    # create output file
    dst_file = dst_dir + grd.name + '_' + dst_varname + '_runoff_bio_NEWS.nc'
    fid_o = nc.Dataset(dst_file, 'w', format='NETCDF3_64BIT')
    fid_o.Description = 'ROMS file'
    fid_o.Author = 'pyroms_toolbox.nc_create_roms_file'
    fid_o.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    fid_o.title = 'ROMS file'

    fid_o.createDimension('xi_rho', np.size(grd.hgrid.mask_rho,1))
    fid_o.createDimension('xi_u', np.size(grd.hgrid.mask_u,1))
    fid_o.createDimension('xi_v', np.size(grd.hgrid.mask_v,1))
    fid_o.createDimension('xi_psi', np.size(grd.hgrid.mask_psi,1))
    fid_o.createDimension('eta_rho', np.size(grd.hgrid.mask_rho,0))
    fid_o.createDimension('eta_u', np.size(grd.hgrid.mask_u,0))
    fid_o.createDimension('eta_v', np.size(grd.hgrid.mask_v,0))
    fid_o.createDimension('eta_psi', np.size(grd.hgrid.mask_psi,0))
    fid_o.createDimension('s_rho', grd.vgrid.N)
    fid_o.createDimension('s_w', grd.vgrid.Np)
    fid_o.createDimension('runoff_time', None)
    fid_o.createVariable('runoff_time', 'f8', ('runoff_time'))
    fid_o.variables['runoff_time'].long_name = 'runoff_time'
    fid_o.variables['runoff_time'].units = 'days since 1900-01-01 00:00:00'
    fid_o.variables['runoff_time'][:] = nctime

    spval = 1.0e+15
    dimensions = ('runoff_time', 'eta_rho', 'xi_rho')
    long_name = tracer_longname[ktr]
    field = dst_varname + ', scalar, series'
    units = tracer_units[ktr]
    fid_o.createVariable(dst_varname, 'f8', dimensions, fill_value=spval)
    fid_o.variables[dst_varname].long_name = long_name
    fid_o.variables[dst_varname].units = units
    fid_o.variables[dst_varname].field = field
    conc = fid_n.variables[list_tracer[ktr]][:]
    for kt in np.arange(nt):
        print('working on timestep', kt)
        runoff = fid_r.variables['Runoff'][kt,:,:]
        river_input = conc * runoff / 1000.
        fid_o.variables[dst_varname][kt,:,:] = river_input


    fid_o.close()

fid_r.close()
fid_n.close()
