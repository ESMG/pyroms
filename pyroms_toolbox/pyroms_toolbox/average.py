# encoding: utf-8

import _average
import types
import pyroms
import numpy as np
import netCDF4
from datetime import datetime

# defined an avg_obj object
class avg_obj(object):
    pass

avg = avg_obj()

def average(var, ncfiles, trange=None, avgfile=None, spval=1e37, timevar='ocean_time'):
    """

    This function computes the temporal average of a given variable or set of
    variables. It uses a shared library object created using f2py to
    calculate this average. Creates a netCDF file if avgfile does not
    equal None, or returns an object containing all of the averages if
    avgfile is equal to None.

    List of Arguments:

    var        - Either a string or list of strings relating to variable
                 names in the netCDF files

    ncfiles    - This variable can be one of three types: a string containing
                 the path to a netCDF file, a string with a Unix wildcard indicator,
                 or a list of strings containing netCDF files.

                 E.g. 1) ncfiles = '/path/to/ncfile.nc'
                      2) ncfiles = '/path/to/nc*.nc'
                      3) ncfiles = ['/path/to/ncfile1.nc','/path/to/ncfile2.nc',...]

    trange      - time range index for ocean time dimension, used if you want to compute
                 average over a given time range.

                 E.g. trange = (10,40)

    avgfile    - A string used for naming a new netCDF file containing only the
                 averages from the variables after being computed.

    timevar    - time variable name. Default is ROMS time variable name "ocean_time".
    """

    if type(var).__name__ == 'list':
        nvar = len(var)
    elif type(var).__name__ == 'str':
        var = [var]
        nvar = len(var)
    else:
        raise ValueError('var must be a str or a list of str')

    avg.ncfiles = pyroms.io.MFDataset(ncfiles)

    ocean_time = pyroms.utility.get_nc_var(timevar, avg.ncfiles)
    Nt = len(ocean_time[:])   

    if trange is None:
        start = 0
        end = Nt
    else:
        if trange[0] >= 0 and trange[1] <= Nt:
            start = trange[0]
            end = min(trange[1]+1, Nt)
        else:
            raise ValueError('trange must be within interval [0, %s].' %Nt)

    print(list(range(start,end)))

    for varname in var:
        name = varname
        vble = pyroms.utility.get_nc_var(varname,avg.ncfiles)
        vsh = vble.shape
        leng = len(vsh)

        # if variable is 4D, enters this conditional
        if leng is 4:
            # create an empty numpy array with dimensions equal to shape of the
            # shape of the variable minus the ocean_time dimension
            incavg = np.zeros((vsh[1],vsh[2],vsh[3]))

            ii = 0
            for i in range(start,end):
                ii = ii + 1
                #calls Fortran function avg3d to perform an incremental average
                incavg = _average.avg3d(vble[i,:],incavg,ii,spval)
            # mask
            incavg = np.ma.masked_values(incavg, spval)
            #sets attribute of avg object to the final temporal average
            setattr(avg, varname, incavg[:])

        # if variable is 3D, enters this conditional
        elif leng is 3:
            # create an empty numpy array with dimensions equal to shape of the
            # shape of the variable minus the ocean_time dimension
            incavg = np.zeros((vsh[1],vsh[2]))

            ii = 0
            for i in range(start,end):
                ii = ii + 1
                #calls Fortran function avg2d to perform an incremental average
                incavg = _average.avg2d(vble[i,:],incavg,ii,spval)
            # mask
            incavg = np.ma.masked_values(incavg, spval)
            #sets attribute of avg object to the final temporal average
            setattr(avg, varname, incavg[:])

        else:
            raise ValueError('Variable must be 3D (time + 2 spacial dims) or \
4D (time + 3 spacial dims)')


    # if avgfile is defined, enter this conditional and begin creating a new netCDF file
    if avgfile is not None:
        # if avgfile is a string, enter this conditional
        if type(avgfile).__name__ == 'str':
            # if original name of avgfile left off .nc in the filename, append it to the string
            if avgfile.find('.nc') == -1:
                avgfile = avgfile+'.nc'

            # creates netCDF file
            nc = netCDF4.Dataset(avgfile, 'w', format='NETCDF3_CLASSIC')
            nc.Description = 'Variable Average File'
            nc.Author = 'pyroms_toolbox.average'
            nc.Created = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            nc.Files = ", ".join(avg.ncfiles._files)

            print('Writing '+str(avgfile)+'...')

            # for each dimension in original netCDF files, recreate the dimensions in the new netCDF file
            for newdim in list(avg.ncfiles.dimensions.keys()):
                if avg.ncfiles.dimensions[newdim].isunlimited():
                    nc.createDimension(newdim,None)
                else:
                    nc.createDimension(newdim,len(avg.ncfiles.dimensions[newdim]))

            tdims = avg.ncfiles.variables[timevar].dimensions
            nc.createVariable(timevar, 'f8', (tdims))
            nc.variables[timevar].long_name = 'averaged time since initialization'
            try:
                nc.variables[timevar].units = avg.ncfiles.variables[timevar].units
                nc.variables[timevar].field = avg.ncfiles.variables[timevar].field
            except:
                nc.variables[timevar].units = 'N/A'
                nc.variables[timevar].field = 'N/A'
            nc.variables[timevar][0] =  ocean_time[start:end].mean()

            # for each variable in var, create a new variable with all dimensions associated with that
            # variable except ocean_time
            for varname in var:
                print('  writting %s...' %varname)

                vardims = avg.ncfiles.variables[varname].dimensions

                nc.createVariable(varname, 'f8', (vardims), fill_value=str(spval))
                nc.variables[varname].long_name = 'Temporal average of variable '+str(varname)

                # try to pull the units from the original netCDF variable
                try:
                    nc.variables[varname].units = avg.ncfiles.variables[varname].units
                    nc.variables[varname].field = avg.ncfiles.variables[varname].field
                # mark as N/A if no units were given to this particular variable
                except:
                    nc.variables[varname].units = 'N/A'
                    nc.variables[varname].field = 'N/A'

                nc.variables[varname][0] = avg.__getattribute__(varname)
            avg.ncfiles.close()
            nc.close()
        else:
            print("avgfile must be a string that equates to the path where this netCDF file is to be placed.")
        return
    else:
        avg.ncfiles.close()
        print("Returning average object...")
        return avg
