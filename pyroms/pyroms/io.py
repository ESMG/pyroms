# encoding: utf-8
'''A thin wrapper for netCDF4.Dataset and netCDF4.MFDataset

This module provides two functions, Dataset and MFDataset, that are similar to the
netCDF[3/4] functions of the same name. This package is a thin wrapper around these
functions, and provides two services. First of all, it will use either netCDF3 or
netCDF4 (prefering the later), so that the netCDF package does not need to be changed
on different systems that only have one or the other. Second, it will pass through
netCDF[3/4] objects unchanged, so that netCDF objects, filenames, lists of files, or
strings with wildcards can be passed to the function indescriminately.

Examples of usage
-----------------

with an input of a string:
    # returns netCDF4.Dataset object based on file
    nc = pyroms.io.Dataset(file)

    # returns MFnetCDF4.Dataset object based on file (with wildcard chars)
    nc = pyroms.io.MFDataset(file)

with an input of a list of files:
    # returns MFDataset object based on list of files
    nc = pyroms.io.Dataset(files)

    # returns MFDataset object based on list of files
    nc = pyroms.io.MFDataset(files)

with an input of a netCDF4.Dataset or MFnetCDF4.Dataset object:
    # passes through netCDF4.Dataset or MFnetCDF4.Dataset object
    nc = pyroms.io.Dataset(nc)

    # passes through MFDataset object based on file (with wildcard chars)
    nc = pyroms.io.MFDataset(nc)
'''
__docformat__ = "restructuredtext en"

from glob import glob

try:
    try:
        import netCDF4 as netCDF
    except:
        import netCDF3 as netCDF

    def Dataset(ncfile):
        """Return an appropriate netcdf object:
                netCDF4 object given a file string
                MFnetCDF4 object given a list of files

            A netCDF4 or MFnetCDF4 object returns itself."""
        if isinstance(ncfile, str):
            return netCDF.Dataset(ncfile, 'r')
        elif isinstance(ncfile, list) or isinstance(ncfile, tuple):
            return netCDF.MFDataset(sorted(ncfile))
        elif hasattr(ncfile, 'variables'):  # accept any oject with a variables attribute
            assert isinstance(ncfile.variables, dict), \
                   'variables attribute must be a dictionary'
            return ncfile
        else:
            raise TypeError('type %s not supported' % type(ncfile))

    Dataset.__doc__ = __doc__

    def MFDataset(ncfile):
        """Return an MFnetCDF4 object given a string or list.  A string is expanded
           with wildcards using glob.  A netCDF4 or MFnetCDF4 object returns itself."""
        if isinstance(ncfile, str):
            ncfiles = glob(ncfile)
            return netCDF.MFDataset(sorted(ncfiles))
        elif isinstance(ncfile, list) or isinstance(ncfile, tuple):
            return netCDF.MFDataset(sorted(ncfile))
        elif hasattr(ncfile, 'variables'):  # accept any oject with a variables attribute
            assert isinstance(ncfile.variables, dict), \
                   'variables attribute must be a dictionary'
            return ncfile
        else:
            raise TypeError('type %s not supported' % type(ncfile))
            return MFnetCDF4.Dataset(files)

    MFDataset.__doc__ = __doc__

except:
    import pyroms.extern.pupynere
    import warnings

    warnings.warn('netCDF[3/4] not found -- using pupynere.')

    def Dataset(ncfile):
        if isinstance(ncfile, str):
            return pupynere.NetCDFFile(ncfile)
        elif isinstance(ncfile, pupynere.NetCDFFile):
            return ncfile
        else:
            raise TypeError('type %s not supported' % type(ncfile))

    Dataset.__doc__ = __doc__


if __name__ == '__main__':
    pass


