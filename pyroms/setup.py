#!/usr/bin/env python

"""pyroms is a suite of tools for working with ROMS.

Requires:
    NumPy (http://numpy.scipy.org)
    matplotlib with the Basemap toolkit (http://matplotlib.sourceforge.net)
    netCDF4 (http://code.google.com/p/netcdf4-python/)

Contains:
    hgrid  -  Tools for dealing with curvilinear grids
      BoundaryInteractor
      _Focus_x
      _Focus_y
       Focus
      CGrid
      CGrid_geo
      Gridgen
      edit_mask_mesh
      get_position_from_map

    vgrid  -  Various vertical coordinates
      s_coordinate
      z_r
      z_w
      z_coordinate

    grid  -  ROMS Grid object
      ROMS_Grid
      ROMS_gridinfo
      print_ROMS_gridinfo
      list_ROMS_gridid
      get_ROMS_hgrid
      get_ROMS_vgrid
      get_ROMS_grid
      write_ROMS_grid                         

    io  -  wrapper for netCDF4
      Dataset
      MFDataset

    cf  -  CF compliant files tools
      time 

    tools  -  Tools specific to the Regional Ocean Modeling System
      roms2z
      z2roms
      zslice
      sslice
      islice
      jslice
      isoslice
      section
      lonslice
      latslice
      section_transport

    utility  -  Some basic tools 
      get_lonlat
      get_ij
      roms_varlist
      get_roms_var
      get_bottom
      get_surface
"""

classifiers = """\
Development Status :: beta
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: MIT
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Software Development :: Libraries :: Python Modules
"""

from numpy.distutils.core import Extension

iso = Extension(name = '_iso',
                sources = ['pyroms/src/iso.f'])

interp = Extension(name = '_interp',
                sources = ['pyroms/src/interp.f'])

obs_interp = Extension(name = '_obs_interp',
                sources = ['pyroms/src/obs_interp.f'])

remapping = Extension(name = '_remapping',
                sources = ['pyroms/src/remapping.f'])

remapping_fast = Extension(name = '_remapping_fast',
                sources = ['pyroms/src/remapping_fast.f'])

remapping_fast_weighted = Extension(name = '_remapping_fast_weighted',
                sources = ['pyroms/src/remapping_fast_weighted.f'])

doclines = __doc__.split("\n")

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(name = "pyroms",
          version = '0.1.0',
          description = doclines[0],
          long_description = "\n".join(doclines[2:]),
          author = "Frederic Castruccio",
          author_email = "frederic@marine.rutgers.edu",
          url = "ftp://marine.rutgers.edu/pub/frederic/pyroms/",
          packages = ['pyroms',
                      'pyroms.remapping',
                      'pyroms.extern'],
          license = 'BSD',
          platforms = ["any"],
          ext_modules = [iso, interp, obs_interp, remapping, remapping_fast, remapping_fast_weighted],
          classifiers = [_f for _f in classifiers.split("\n") if _f],
          )
